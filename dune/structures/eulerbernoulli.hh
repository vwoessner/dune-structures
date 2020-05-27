#ifndef DUNE_STRUCTURES_EULERBERNOULLI_HH
#define DUNE_STRUCTURES_EULERBERNOULLI_HH

/** Local operators to implement the beam-truss elements for
 *  Euler-Bernoulli beams used in
 *
 *  Hansbo P., Larson M. and Larsson K.
 *  Cut Finite Element Methods for Linear Elasticity Problems
 */

#include<dune/geometry/affinegeometry.hh>
#include<dune/pdelab.hh>
#include<dune/structures/elasticity.hh>
#include<dune/structures/material.hh>
#include<array>
#include<map>
#include<tuple>


template<int dim>
class FibreParametrizationBase
{
  public:
  using GlobalCoordinate = Dune::FieldVector<double, dim>;
  virtual ~FibreParametrizationBase() = default;
  virtual GlobalCoordinate eval(double t) const = 0;
  virtual GlobalCoordinate tangent(double t) const = 0;
};


template<int dim>
class StraightFibre
  : public FibreParametrizationBase<dim>
{
  public:
  using GlobalCoordinate = Dune::FieldVector<double, dim>;
  virtual ~StraightFibre() = default;

  StraightFibre(const Dune::ParameterTree& param)
    : start(param.get<GlobalCoordinate>("start"))
    , dir(param.get<GlobalCoordinate>("end") - start)
  {}

  virtual GlobalCoordinate eval(double t) const override
  {
    return start + t * dir;
  }

  virtual GlobalCoordinate tangent(double t) const override
  {
    GlobalCoordinate tang(dir);
    tang /= tang.two_norm();
    return tang;
  }

  private:
  GlobalCoordinate start, dir;
};


template<typename LocalBasis>
class BasisEvaluator
{
  public:
  using RFT = typename LocalBasis::Traits::RangeFieldType;

  BasisEvaluator(const LocalBasis& basis) : basis(basis)
  {}

  void update(const typename LocalBasis::Traits::DomainType& coord)
  {
    basis.partial({0, 0}, coord, phi);
    basis.partial({1, 0}, coord, d1phi[0]);
    basis.partial({0, 1}, coord, d1phi[1]);
    basis.partial({2, 0}, coord, d2phi[0]);
    basis.partial({1, 1}, coord, d2phi[1]);
    basis.partial({0, 2}, coord, d2phi[2]);
  }

  RFT function(std::size_t i) const
  {
    return phi[i][0];
  }

  RFT jacobian(std::size_t i, std::size_t j) const
  {
    return d1phi[j][i][0];
  }

  RFT hessian(std::size_t i, std::size_t j, std::size_t k) const
  {
    return d2phi[j + k][i][0];
  }

  private:
  const LocalBasis& basis;
  std::vector<typename LocalBasis::Traits::RangeType> phi;
  std::array<std::vector<typename LocalBasis::Traits::RangeType>, 2> d1phi;
  std::array<std::vector<typename LocalBasis::Traits::RangeType>, 3> d2phi;
};


/** An operator that implement a single Euler-Bernoulli beam
 *
 * For multiple beams, use a sum of operators. This may not be super
 * high performance, but things are already complicated enough!
 */
template<typename GFS>
class EulerBernoulli2DLocalOperator
  :  public Dune::PDELab::FullSkeletonPattern
  ,  public Dune::PDELab::FullVolumePattern
  ,  public Dune::PDELab::LocalOperatorDefaultFlags
  ,  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  ,  public Dune::PDELab::NumericalJacobianVolume<EulerBernoulli2DLocalOperator<GFS>>
  ,  public Dune::PDELab::NumericalJacobianSkeleton<EulerBernoulli2DLocalOperator<GFS>>
{
  public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  enum { doPatternSkeleton = true };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doAlphaSkeleton  = true };

  EulerBernoulli2DLocalOperator(const Dune::ParameterTree& rootparams, const Dune::ParameterTree& params, std::shared_ptr<const GFS> gfs)
    : Dune::PDELab::NumericalJacobianVolume<EulerBernoulli2DLocalOperator<GFS>>(1e-9)
    , Dune::PDELab::NumericalJacobianSkeleton<EulerBernoulli2DLocalOperator<GFS>>(1e-9)
    , is(gfs->gridView().indexSet())
  {
    // Some debugging switches to reduce recompilations in debugging
    bool verbose = params.get<bool>("verbose", false);

	// Extract stabilization parameter
	beta = params.get<double>("stabilization_parameter", 1.0);

    // Parse fibres from the configuration
    auto fibrestr = rootparams.get<std::string>("grid.fibres.fibres", "");
    auto fibres = str_split(fibrestr);
    for (auto fibre: fibres)
    {
      str_trim(fibre);
      auto fibreconfig = rootparams.sub("grid.fibres").sub(fibre);

      if (fibreconfig.get<std::string>("shape") == "cylinder")
        fibre_parametrizations.push_back(std::make_shared<StraightFibre<2>>(fibreconfig));
      else
        DUNE_THROW(Dune::Exception, "Fibre shape not supported!");

      // Parse the Youngs modulus of the fibre.
      auto modulus = fibreconfig.get<double>("youngs_modulus", 0.0);
      if (modulus == 0.0)
        std::cout << "youngs_modulus not set for all fibres!" << std::endl;
      fibre_modulus.push_back(modulus);

      // Parse the fibre radii
      auto radius = fibreconfig.get<double>("radius", 0.0);
      if (radius == 0.0)
        std::cout << "radius not set for all fibres!" << std::endl;
      fibre_radii.push_back(radius);
    }
    std::cout << "Parsed a total of " << fibre_parametrizations.size() << " fibres from the configuration." << std::endl;

    /* Find intersections of fibres with the grid. The outline is rougly the following:
       - Find the first simplex that the fibre intersects by doing hierarchical search on fibre(0)
       - If fibre(0) is outside the grid, do a more expensive method to find the first simplex
       - Do a bisection of the fibre interval to find the parameter t where it leaves the simplex
       - Store the parameter interval and the simplex index
       - Go to the neighboring simplex and repeat until at a boundary intersection or at fibre(1).
    */
    for (std::size_t fibindex = 0; fibindex < fibre_parametrizations.size(); ++fibindex)
    {
      auto fibre = fibre_parametrizations[fibindex];
      auto gv = gfs->gridView();
      int index = -1;
      typename GFS::Traits::GridView::template Codim<0>::Entity element = *(gv.template begin<0>());

      // This tolerance threshold is used in the following
      double tol = 1e-12;

      // We iterate this algorithm, because in rare scenarios the strategy of going from cell to cell
      // across intersections fails. This is when the curve traverses through a grid vertex from one cell
      // to a non-adjacent one. In these cases we restart the algorithm by doing a grid search for the
      // current cell.
      double t = 0.0;
      while(t + tol < 1.0)
      {
        if (t > tol)
          std::cout << "A restart of the curve cutting algorithm was needed at t=" << t << std::endl;

        // Find the first cell in the grid that the curve traverses
        for(const auto& e : elements(gv))
        {
          auto geo = e.geometry();
          if (referenceElement(geo).checkInside(geo.local(fibre->eval(t + tol))))
          {
            index = is.index(e);
            element = e;
            break;
          }
        }

        // Coordinate not in the grid, but this is not the start => It is the end.
        if ((index == -1) && (t > tol))
          break;

        // The coordinate was not found in the grid, we need to bisect the curve
        // to find a starting element within the grid.
        if (index == -1)
        {
          double bisect_a = t;
          double bisect_b = 0.5;
          while (bisect_b - bisect_a > tol)
          {
            double mid = 0.5 * (bisect_a + bisect_b);
            index = -1;
            for(const auto& e : elements(gv))
            {
              auto geo = e.geometry();
              if (referenceElement(geo).checkInside(geo.local(fibre->eval(mid))))
              {
                index = is.index(e);
                element = e;
                break;
              }
            }

            if (index == -1)
              bisect_a = mid;
            else
              bisect_b = mid;
          }

          // Restart the outer loop
          t = bisect_a;
          continue;
        }

        // Go from cell to cell by finding the intersection that the
        while(index != -1)
        {
          // Find the value of t where the curve leaves the element
          double bisect_a = t;
          double bisect_b = 1.0;
          while (bisect_b - bisect_a > tol)
          {
            double mid = 0.5 * (bisect_a + bisect_b);
            auto geo = element.geometry();
            if (referenceElement(geo).checkInside(geo.local(fibre->eval(mid))))
              bisect_a = mid;
            else
              bisect_b = mid;
          }

          // Store the fibre segment that we found
          element_fibre_intersections[index].push_back({fibindex, t, bisect_a});
          t = bisect_a;

          // Find the intersection and the outside cell
          index = -1;
          for (auto intersection : intersections(gv, element))
          {
            if (!intersection.boundary())
            {
              auto geo = intersection.outside().geometry();
              if (referenceElement(geo).checkInside(geo.local(fibre->eval(t + tol))))
              {
                int inside_index = is.index(element);
                element = intersection.outside();
                index = is.index(element);
                face_fibre_intersections[std::make_pair(inside_index, index)].push_back({fibindex, t});
                break;
              }
            }
          }
        }
      }
    }

    if (verbose)
    {
      std::cout << "Fibre intersection summary:" << std::endl;
      for (auto [cell, info] : element_fibre_intersections)
        for (auto sinfo : info)
          std::cout << "Cell " << cell << " intersects fibre " << std::get<0>(sinfo)
                    << " on the curve interval [" << std::get<1>(sinfo) << "," << std::get<2>(sinfo) << "]" << std::endl;
      for (auto [cells, info] : face_fibre_intersections)
      {
        auto [inside, outside] = cells;
        for (auto sinfo : info)
          std::cout << "Facet between cell " << inside << " and " << outside << " intersects fibre "
                    << std::get<0>(sinfo) << " at t=" << std::get<1>(sinfo) << std::endl;
      }
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // See whether something needs to be done on this cell
    auto entity = eg.entity();
    auto it = element_fibre_intersections.find(is.index(entity));
    if (it == element_fibre_intersections.end())
      return;

    // Extract some necessary information
    using namespace Dune::Indices;
    auto child_0 = child(lfsu, _0);
    auto child_1 = child(lfsu, _1);
    auto cellgeo = entity.geometry();

    BasisEvaluator<typename LFSU::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType> basis(child_0.finiteElement().localBasis());

    // Get the force coefficients
    coefficient_force_lfs->bind(eg.entity());
    typename CoefficientForceVector::template LocalView<CoefficientForceLFSCache> coefficient_force_view(*coefficient_force_vector);
    coefficient_force_lfs_cache->update();
    Dune::PDELab::LocalVector<double> local_coefficient_force_vector(coefficient_force_lfs->size());
    coefficient_force_view.bind(*coefficient_force_lfs_cache);
    coefficient_force_view.read(local_coefficient_force_vector);
    coefficient_force_view.unbind();

    // Loop over curve segments in this cell. most of the time this will be just one.
    auto segments = it->second;
    for (auto segment : segments)
    {
      auto fibre = fibre_parametrizations[std::get<0>(segment)];
      auto start = fibre->eval(std::get<1>(segment));
      auto stop = fibre->eval(std::get<2>(segment));

      // Construct the geometry of the 1D inclusion embedded into the reference element of the cell
      using LineGeometry = Dune::AffineGeometry<double, 1, 2>;
      const auto& linerefelem = Dune::Geo::ReferenceElements<double, 1>::simplex();
      LineGeometry linegeo(linerefelem, std::vector<Dune::FieldVector<double, 2>>{cellgeo.local(start), cellgeo.local(stop)});
      LineGeometry global_linegeo(linerefelem, std::vector<Dune::FieldVector<double, 2>>{start, stop});

      // Iterate over the quadrature points - currently always midpoint rule
      for (const auto& ip : quadratureRule(linegeo, 2))
      {
        // Position in reference element of the cell
        auto pos = linegeo.global(ip.position());

        // Evaluate the basis, its jacobian and its hessians
        basis.update(pos);

        // Evaluate the displacement field u
        Dune::FieldVector<double, 2> u(0.0);
        for (std::size_t c=0; c<2; ++c)
          for (std::size_t i=0; i<child_0.size(); i++)
            u[c] += x(child(lfsu, c), i) * basis.function(i);

        // And the jacobian of displacement
        Dune::FieldMatrix<double, 2, 2> d1u(0.0);
        for (std::size_t c=0; c<2; ++c)
          for (std::size_t d=0; d<2; ++d)
            for (std::size_t i=0; i<child_0.size(); ++i)
              d1u[c][d] += x(child(lfsu, c), i) * basis.jacobian(i, d);

        // And finally the hessian of the displacement
        std::array<Dune::FieldMatrix<double, 2, 2>, 2> d2u{0.0, 0.0};
        for (std::size_t c=0; c<2; ++c)
          for (std::size_t d0=0; d0<2; ++d0)
            for (std::size_t d1=0; d1<2; ++d1)
              for (std::size_t i=0; i<child_0.size(); ++i)
                d2u[c][d0][d1] += x(child(lfsu, c), i) * basis.hessian(i, d0, d1);

        // Evaluate the jacobian inverse transposed
        auto jit = cellgeo.jacobianInverseTransposed(pos);

        // The tangential vector for the curve
        // Determination of the evaluation parameter here is a bit flaky.
        auto t = fibre->tangent(std::get<1>(segment) + ip.position() * (std::get<2>(segment) - std::get<1>(segment)));

        // Evaluate the body force vector
        Dune::FieldVector<double, 2> force(0.0);
        for (std::size_t k = 0; k < 2; ++k)
          for (std::size_t i = 0; i < child_0.size(); ++i)
            force[k] += local_coefficient_force_vector(lfsv.child(k), i) * basis.function(i);

        // Extract physical parameter of the fibre
        auto E = fibre_modulus[std::get<0>(segment)];
        auto d = fibre_radii[std::get<0>(segment)];
        auto A = d;
        auto I = (d*d*d) / 12.0;

        // The needed tangential derivative quantities. These expressions are generated
        // using the generate_tangential_derivatives Python script.
        auto dtut = ((d1u[1][1] * jit[1][1] + d1u[1][0] * jit[1][0]) * t[1] + (d1u[0][1] * jit[1][1] + d1u[0][0] * jit[1][0]) * t[0]) * t[1] + ((d1u[1][1] * jit[0][1] + d1u[1][0] * jit[0][0]) * t[1] + (d1u[0][1] * jit[0][1] + d1u[0][0] * jit[0][0]) * t[0]) * t[0];
        auto dt2un = ((t[0] * (jit[1][1] * (jit[1][1] * d2u[1][1][1] + jit[1][0] * d2u[1][1][0]) + jit[1][0] * (jit[1][1] * d2u[1][0][1] + jit[1][0] * d2u[1][0][0])) + (jit[1][1] * (jit[1][1] * d2u[0][1][1] + jit[1][0] * d2u[0][1][0]) + jit[1][0] * (jit[1][1] * d2u[0][0][1] + jit[1][0] * d2u[0][0][0])) * (-1) * t[1]) * t[1] + (t[0] * (jit[0][1] * (jit[1][1] * d2u[1][1][1] + jit[1][0] * d2u[1][1][0]) + jit[0][0] * (jit[1][1] * d2u[1][0][1] + jit[1][0] * d2u[1][0][0])) + (jit[0][1] * (jit[1][1] * d2u[0][1][1] + jit[1][0] * d2u[0][1][0]) + jit[0][0] * (jit[1][1] * d2u[0][0][1] + jit[1][0] * d2u[0][0][0])) * (-1) * t[1]) * t[0]) * t[1] + ((t[0] * (jit[1][1] * (jit[0][1] * d2u[1][1][1] + jit[0][0] * d2u[1][1][0]) + jit[1][0] * (jit[0][1] * d2u[1][0][1] + jit[0][0] * d2u[1][0][0])) + (jit[1][1] * (jit[0][1] * d2u[0][1][1] + jit[0][0] * d2u[0][1][0]) + jit[1][0] * (jit[0][1] * d2u[0][0][1] + jit[0][0] * d2u[0][0][0])) * (-1) * t[1]) * t[1] + (t[0] * (jit[0][1] * (jit[0][1] * d2u[1][1][1] + jit[0][0] * d2u[1][1][0]) + jit[0][0] * (jit[0][1] * d2u[1][0][1] + jit[0][0] * d2u[1][0][0])) + (jit[0][1] * (jit[0][1] * d2u[0][1][1] + jit[0][0] * d2u[0][1][0]) + jit[0][0] * (jit[0][1] * d2u[0][0][1] + jit[0][0] * d2u[0][0][0])) * (-1) * t[1]) * t[0]) * t[0];

        for (std::size_t i=0; i<child_0.size(); ++i)
        {
          auto dtvt_0 = t[1] * (jit[1][1] * basis.jacobian(i, 1) + jit[1][0] * basis.jacobian(i, 0)) * t[0] + t[0] * (jit[0][1] * basis.jacobian(i, 1) + jit[0][0] * basis.jacobian(i, 0)) * t[0];
          auto dtvt_1 = t[1] * (jit[1][1] * basis.jacobian(i, 1) + jit[1][0] * basis.jacobian(i, 0)) * t[1] + t[0] * (jit[0][1] * basis.jacobian(i, 1) + jit[0][0] * basis.jacobian(i, 0)) * t[1];
          auto dt2vn_0 = (t[1] * (jit[1][1] * (jit[1][1] * basis.hessian(i, 1, 1) + jit[1][0] * basis.hessian(i, 1, 0)) + jit[1][0] * (jit[1][1] * basis.hessian(i, 0, 1) + jit[1][0] * basis.hessian(i, 0, 0))) * (-1) * t[1] + t[0] * (jit[0][1] * (jit[1][1] * basis.hessian(i, 1, 1) + jit[1][0] * basis.hessian(i, 1, 0)) + jit[0][0] * (jit[1][1] * basis.hessian(i, 0, 1) + jit[1][0] * basis.hessian(i, 0, 0))) * (-1) * t[1]) * t[1] + (t[1] * (jit[1][1] * (jit[0][1] * basis.hessian(i, 1, 1) + jit[0][0] * basis.hessian(i, 1, 0)) + jit[1][0] * (jit[0][1] * basis.hessian(i, 0, 1) + jit[0][0] * basis.hessian(i, 0, 0))) * (-1) * t[1] + t[0] * (jit[0][1] * (jit[0][1] * basis.hessian(i, 1, 1) + jit[0][0] * basis.hessian(i, 1, 0)) + jit[0][0] * (jit[0][1] * basis.hessian(i, 0, 1) + jit[0][0] * basis.hessian(i, 0, 0))) * (-1) * t[1]) * t[0];
          auto dt2vn_1 = (t[1] * t[0] * (jit[1][1] * (jit[1][1] * basis.hessian(i, 1, 1) + jit[1][0] * basis.hessian(i, 1, 0)) + jit[1][0] * (jit[1][1] * basis.hessian(i, 0, 1) + jit[1][0] * basis.hessian(i, 0, 0))) + t[0] * t[0] * (jit[0][1] * (jit[1][1] * basis.hessian(i, 1, 1) + jit[1][0] * basis.hessian(i, 1, 0)) + jit[0][0] * (jit[1][1] * basis.hessian(i, 0, 1) + jit[1][0] * basis.hessian(i, 0, 0)))) * t[1] + (t[1] * t[0] * (jit[1][1] * (jit[0][1] * basis.hessian(i, 1, 1) + jit[0][0] * basis.hessian(i, 1, 0)) + jit[1][0] * (jit[0][1] * basis.hessian(i, 0, 1) + jit[0][0] * basis.hessian(i, 0, 0))) + t[0] * t[0] * (jit[0][1] * (jit[0][1] * basis.hessian(i, 1, 1) + jit[0][0] * basis.hessian(i, 1, 0)) + jit[0][0] * (jit[0][1] * basis.hessian(i, 0, 1) + jit[0][0] * basis.hessian(i, 0, 0)))) * t[0];

          r.accumulate(lfsu.child(0), i, (E*A*dtut*dtvt_0 + E*I*dt2un*dt2vn_0 - A*force[0]*basis.function(i)) * ip.weight() * global_linegeo.integrationElement(ip.position()));
          r.accumulate(lfsu.child(1), i, (E*A*dtut*dtvt_1 + E*I*dt2un*dt2vn_1 - A*force[1]*basis.function(i)) * ip.weight() * global_linegeo.integrationElement(ip.position()));
        }
      }
    }
  }

  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton (const IG& ig,
                       const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                       const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                       R& r_s, R& r_n) const
  {
    // The notion of inside and outside is quite complex in this situation.
    // We are traversing the 2D grid, whose intersections have inside and outside.
    // But these do not necessarily match the notion of inside and outside of
    // the 1D discretization. We therefore need to remap these references in some
    // cases.
    bool flipped = false;

    auto it = face_fibre_intersections.find(std::make_pair(is.index(ig.inside()), is.index(ig.outside())));
    if (it == face_fibre_intersections.end())
    {
      it = face_fibre_intersections.find(std::make_pair(is.index(ig.outside()), is.index(ig.inside())));
      if (it == face_fibre_intersections.end())
        return;

      flipped = true;
    }

    // Extract some necessary information
    using namespace Dune::Indices;
    auto child_0 = child(lfsu_s, _0);
    auto child_1 = child(lfsu_s, _1);
    auto cellgeo_s = (flipped ? ig.outside() : ig.inside()).geometry();
    auto cellgeo_n = (flipped ? ig.inside() : ig.outside()).geometry();

    using Evaluator = BasisEvaluator<typename LFSU::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType>;
    Evaluator basis_s(child_0.finiteElement().localBasis());
    Evaluator basis_n(child_0.finiteElement().localBasis());

    // Loop over curve segments in this cell. most of the time this will be just one.
    auto segments = it->second;
    for (auto segment : segments)
    {
      auto fibre = fibre_parametrizations[std::get<0>(segment)];
      auto tparam = fibre->eval(std::get<1>(segment));

      // Position in reference element of the inside/outside cell
      auto pos_s = cellgeo_s.local(tparam);
      auto pos_n = cellgeo_n.local(tparam);

      // Evaluate the basis, its jacobian and its hessians in inside/outside cell
      basis_s.update(pos_s);
      basis_n.update(pos_n);

      // Evaluate the displacement field u
      Dune::FieldVector<double, 2> u_s(0.0), u_n(0.0);
      for (std::size_t c=0; c<2; ++c)
        for (std::size_t i=0; i<child_0.size(); i++)
        {
          u_s[c] += (flipped ? x_n : x_s)(child(lfsu_s, c), i) * basis_s.function(i);
          u_n[c] += (flipped ? x_s : x_n)(child(lfsu_n, c), i) * basis_n.function(i);
        }

      // And the jacobian of displacement
      Dune::FieldMatrix<double, 2, 2> d1u_s(0.0), d1u_n(0.0);
      for (std::size_t c=0; c<2; ++c)
        for (std::size_t d=0; d<2; ++d)
          for (std::size_t i=0; i<child_0.size(); ++i)
          {
            d1u_s[c][d] += (flipped ? x_n : x_s)(child(lfsu_s, c), i) * basis_s.jacobian(i, d);
            d1u_n[c][d] += (flipped ? x_s : x_n)(child(lfsu_n, c), i) * basis_n.jacobian(i, d);
          }

      // And finally the hessian of the displacement
      std::array<Dune::FieldMatrix<double, 2, 2>, 2> d2u_s{0.0, 0.0}, d2u_n{0.0, 0.0};
      for (std::size_t c=0; c<2; ++c)
        for (std::size_t d0=0; d0<2; ++d0)
          for (std::size_t d1=0; d1<2; ++d1)
            for (std::size_t i=0; i<child_0.size(); ++i)
            {
              d2u_s[c][d0][d1] += (flipped ? x_n : x_s)(child(lfsu_s, c), i) * basis_s.hessian(i, d0, d1);
              d2u_n[c][d0][d1] += (flipped ? x_s : x_n)(child(lfsu_n, c), i) * basis_n.hessian(i, d0, d1);
            }

      // Evaluate the jacobian inverse transposed in inside/outside cell
      auto jit_s = cellgeo_s.jacobianInverseTransposed(pos_s);
      auto jit_n = cellgeo_n.jacobianInverseTransposed(pos_n);

      // The tangential vector for the curve
      auto t = fibre->tangent(std::get<1>(segment));

      // Extract physical parameter of the fibre
      auto E = fibre_modulus[std::get<0>(segment)];
      auto d = fibre_radii[std::get<0>(segment)];
      auto A = d;
      auto I = (d*d*d) / 12.0;

      // The needed tangential derivative quantities. These expressions are generated
      // using the generate_tangential_derivatives Python script.
      auto sk_dtun = (t[0] * (jit_s[1][1] * d1u_s[1][1] + jit_s[1][0] * d1u_s[1][0]) + (jit_s[1][1] * d1u_s[0][1] + jit_s[1][0] * d1u_s[0][0]) * (-1) * t[1]) * t[1] + (t[0] * (jit_s[0][1] * d1u_s[1][1] + jit_s[0][0] * d1u_s[1][0]) + (jit_s[0][1] * d1u_s[0][1] + jit_s[0][0] * d1u_s[0][0]) * (-1) * t[1]) * t[0] - ((t[0] * (jit_n[1][1] * d1u_n[1][1] + jit_n[1][0] * d1u_n[1][0]) + (jit_n[1][1] * d1u_n[0][1] + jit_n[1][0] * d1u_n[0][0]) * (-1) * t[1]) * t[1] + (t[0] * (jit_n[0][1] * d1u_n[1][1] + jit_n[0][0] * d1u_n[1][0]) + (jit_n[0][1] * d1u_n[0][1] + jit_n[0][0] * d1u_n[0][0]) * (-1) * t[1]) * t[0]);
      auto sk_dt2un = 0.5 * (((t[0] * ((jit_s[1][1] * d2u_s[1][1][1] + jit_s[1][0] * d2u_s[1][1][0]) * jit_s[1][1] + (jit_s[1][1] * d2u_s[1][0][1] + jit_s[1][0] * d2u_s[1][0][0]) * jit_s[1][0]) + ((jit_s[1][1] * d2u_s[0][1][1] + jit_s[1][0] * d2u_s[0][1][0]) * jit_s[1][1] + (jit_s[1][1] * d2u_s[0][0][1] + jit_s[1][0] * d2u_s[0][0][0]) * jit_s[1][0]) * (-1) * t[1]) * t[1] + (t[0] * ((jit_s[1][1] * d2u_s[1][1][1] + jit_s[1][0] * d2u_s[1][1][0]) * jit_s[0][1] + (jit_s[1][1] * d2u_s[1][0][1] + jit_s[1][0] * d2u_s[1][0][0]) * jit_s[0][0]) + ((jit_s[1][1] * d2u_s[0][1][1] + jit_s[1][0] * d2u_s[0][1][0]) * jit_s[0][1] + (jit_s[1][1] * d2u_s[0][0][1] + jit_s[1][0] * d2u_s[0][0][0]) * jit_s[0][0]) * (-1) * t[1]) * t[0]) * t[1] + ((t[0] * ((jit_s[0][1] * d2u_s[1][1][1] + jit_s[0][0] * d2u_s[1][1][0]) * jit_s[1][1] + (jit_s[0][1] * d2u_s[1][0][1] + jit_s[0][0] * d2u_s[1][0][0]) * jit_s[1][0]) + ((jit_s[0][1] * d2u_s[0][1][1] + jit_s[0][0] * d2u_s[0][1][0]) * jit_s[1][1] + (jit_s[0][1] * d2u_s[0][0][1] + jit_s[0][0] * d2u_s[0][0][0]) * jit_s[1][0]) * (-1) * t[1]) * t[1] + (t[0] * ((jit_s[0][1] * d2u_s[1][1][1] + jit_s[0][0] * d2u_s[1][1][0]) * jit_s[0][1] + (jit_s[0][1] * d2u_s[1][0][1] + jit_s[0][0] * d2u_s[1][0][0]) * jit_s[0][0]) + ((jit_s[0][1] * d2u_s[0][1][1] + jit_s[0][0] * d2u_s[0][1][0]) * jit_s[0][1] + (jit_s[0][1] * d2u_s[0][0][1] + jit_s[0][0] * d2u_s[0][0][0]) * jit_s[0][0]) * (-1) * t[1]) * t[0]) * t[0] + ((t[0] * ((jit_n[1][1] * d2u_n[1][1][1] + jit_n[1][0] * d2u_n[1][1][0]) * jit_n[1][1] + (jit_n[1][1] * d2u_n[1][0][1] + jit_n[1][0] * d2u_n[1][0][0]) * jit_n[1][0]) + ((jit_n[1][1] * d2u_n[0][1][1] + jit_n[1][0] * d2u_n[0][1][0]) * jit_n[1][1] + (jit_n[1][1] * d2u_n[0][0][1] + jit_n[1][0] * d2u_n[0][0][0]) * jit_n[1][0]) * (-1) * t[1]) * t[1] + (t[0] * ((jit_n[1][1] * d2u_n[1][1][1] + jit_n[1][0] * d2u_n[1][1][0]) * jit_n[0][1] + (jit_n[1][1] * d2u_n[1][0][1] + jit_n[1][0] * d2u_n[1][0][0]) * jit_n[0][0]) + ((jit_n[1][1] * d2u_n[0][1][1] + jit_n[1][0] * d2u_n[0][1][0]) * jit_n[0][1] + (jit_n[1][1] * d2u_n[0][0][1] + jit_n[1][0] * d2u_n[0][0][0]) * jit_n[0][0]) * (-1) * t[1]) * t[0]) * t[1] + ((t[0] * ((jit_n[0][1] * d2u_n[1][1][1] + jit_n[0][0] * d2u_n[1][1][0]) * jit_n[1][1] + (jit_n[0][1] * d2u_n[1][0][1] + jit_n[0][0] * d2u_n[1][0][0]) * jit_n[1][0]) + ((jit_n[0][1] * d2u_n[0][1][1] + jit_n[0][0] * d2u_n[0][1][0]) * jit_n[1][1] + (jit_n[0][1] * d2u_n[0][0][1] + jit_n[0][0] * d2u_n[0][0][0]) * jit_n[1][0]) * (-1) * t[1]) * t[1] + (t[0] * ((jit_n[0][1] * d2u_n[1][1][1] + jit_n[0][0] * d2u_n[1][1][0]) * jit_n[0][1] + (jit_n[0][1] * d2u_n[1][0][1] + jit_n[0][0] * d2u_n[1][0][0]) * jit_n[0][0]) + ((jit_n[0][1] * d2u_n[0][1][1] + jit_n[0][0] * d2u_n[0][1][0]) * jit_n[0][1] + (jit_n[0][1] * d2u_n[0][0][1] + jit_n[0][0] * d2u_n[0][0][0]) * jit_n[0][0]) * (-1) * t[1]) * t[0]) * t[0]);

      for (std::size_t i=0; i<child_0.size(); ++i)
      {
        auto dtvn_n_0 = t[1] * (jit_n[1][1] * basis_n.jacobian(i, 1) + jit_n[1][0] * basis_n.jacobian(i, 0)) * (-1) * t[1] + t[0] * (jit_n[0][1] * basis_n.jacobian(i, 1) + jit_n[0][0] * basis_n.jacobian(i, 0)) * (-1) * t[1];
        auto dtvn_n_1 = t[1] * t[0] * (jit_n[1][1] * basis_n.jacobian(i, 1) + jit_n[1][0] * basis_n.jacobian(i, 0)) + t[0] * t[0] * (jit_n[0][1] * basis_n.jacobian(i, 1) + jit_n[0][0] * basis_n.jacobian(i, 0));

        auto dtvn_s_0 = t[1] * (jit_s[1][1] * basis_s.jacobian(i, 1) + jit_s[1][0] * basis_s.jacobian(i, 0)) * (-1) * t[1] + t[0] * (jit_s[0][1] * basis_s.jacobian(i, 1) + jit_s[0][0] * basis_s.jacobian(i, 0)) * (-1) * t[1];
        auto dtvn_s_1 = t[1] * t[0] * (jit_s[1][1] * basis_s.jacobian(i, 1) + jit_s[1][0] * basis_s.jacobian(i, 0)) + t[0] * t[0] * (jit_s[0][1] * basis_s.jacobian(i, 1) + jit_s[0][0] * basis_s.jacobian(i, 0));

        auto dt2vn_n_0 = (t[1] * (jit_n[1][1] * (jit_n[1][1] * basis_n.hessian(i, 1, 1) + jit_n[1][0] * basis_n.hessian(i, 1, 0)) + jit_n[1][0] * (jit_n[1][1] * basis_n.hessian(i, 0, 1) + jit_n[1][0] * basis_n.hessian(i, 0, 0))) * (-1) * t[1] + t[0] * (jit_n[0][1] * (jit_n[1][1] * basis_n.hessian(i, 1, 1) + jit_n[1][0] * basis_n.hessian(i, 1, 0)) + jit_n[0][0] * (jit_n[1][1] * basis_n.hessian(i, 0, 1) + jit_n[1][0] * basis_n.hessian(i, 0, 0))) * (-1) * t[1]) * t[1] + (t[1] * (jit_n[1][1] * (jit_n[0][1] * basis_n.hessian(i, 1, 1) + jit_n[0][0] * basis_n.hessian(i, 1, 0)) + jit_n[1][0] * (jit_n[0][1] * basis_n.hessian(i, 0, 1) + jit_n[0][0] * basis_n.hessian(i, 0, 0))) * (-1) * t[1] + t[0] * (jit_n[0][1] * (jit_n[0][1] * basis_n.hessian(i, 1, 1) + jit_n[0][0] * basis_n.hessian(i, 1, 0)) + jit_n[0][0] * (jit_n[0][1] * basis_n.hessian(i, 0, 1) + jit_n[0][0] * basis_n.hessian(i, 0, 0))) * (-1) * t[1]) * t[0];
        auto dt2vn_n_1 = (t[1] * t[0] * (jit_n[1][1] * (jit_n[1][1] * basis_n.hessian(i, 1, 1) + jit_n[1][0] * basis_n.hessian(i, 1, 0)) + jit_n[1][0] * (jit_n[1][1] * basis_n.hessian(i, 0, 1) + jit_n[1][0] * basis_n.hessian(i, 0, 0))) + t[0] * t[0] * (jit_n[0][1] * (jit_n[1][1] * basis_n.hessian(i, 1, 1) + jit_n[1][0] * basis_n.hessian(i, 1, 0)) + jit_n[0][0] * (jit_n[1][1] * basis_n.hessian(i, 0, 1) + jit_n[1][0] * basis_n.hessian(i, 0, 0)))) * t[1] + (t[1] * t[0] * (jit_n[1][1] * (jit_n[0][1] * basis_n.hessian(i, 1, 1) + jit_n[0][0] * basis_n.hessian(i, 1, 0)) + jit_n[1][0] * (jit_n[0][1] * basis_n.hessian(i, 0, 1) + jit_n[0][0] * basis_n.hessian(i, 0, 0))) + t[0] * t[0] * (jit_n[0][1] * (jit_n[0][1] * basis_n.hessian(i, 1, 1) + jit_n[0][0] * basis_n.hessian(i, 1, 0)) + jit_n[0][0] * (jit_n[0][1] * basis_n.hessian(i, 0, 1) + jit_n[0][0] * basis_n.hessian(i, 0, 0)))) * t[0];

        auto dt2vn_s_0 = (t[1] * (jit_s[1][1] * (jit_s[1][1] * basis_s.hessian(i, 1, 1) + jit_s[1][0] * basis_s.hessian(i, 1, 0)) + jit_s[1][0] * (jit_s[1][1] * basis_s.hessian(i, 0, 1) + jit_s[1][0] * basis_s.hessian(i, 0, 0))) * (-1) * t[1] + t[0] * (jit_s[0][1] * (jit_s[1][1] * basis_s.hessian(i, 1, 1) + jit_s[1][0] * basis_s.hessian(i, 1, 0)) + jit_s[0][0] * (jit_s[1][1] * basis_s.hessian(i, 0, 1) + jit_s[1][0] * basis_s.hessian(i, 0, 0))) * (-1) * t[1]) * t[1] + (t[1] * (jit_s[1][1] * (jit_s[0][1] * basis_s.hessian(i, 1, 1) + jit_s[0][0] * basis_s.hessian(i, 1, 0)) + jit_s[1][0] * (jit_s[0][1] * basis_s.hessian(i, 0, 1) + jit_s[0][0] * basis_s.hessian(i, 0, 0))) * (-1) * t[1] + t[0] * (jit_s[0][1] * (jit_s[0][1] * basis_s.hessian(i, 1, 1) + jit_s[0][0] * basis_s.hessian(i, 1, 0)) + jit_s[0][0] * (jit_s[0][1] * basis_s.hessian(i, 0, 1) + jit_s[0][0] * basis_s.hessian(i, 0, 0))) * (-1) * t[1]) * t[0];
        auto dt2vn_s_1 = (t[1] * t[0] * (jit_s[1][1] * (jit_s[1][1] * basis_s.hessian(i, 1, 1) + jit_s[1][0] * basis_s.hessian(i, 1, 0)) + jit_s[1][0] * (jit_s[1][1] * basis_s.hessian(i, 0, 1) + jit_s[1][0] * basis_s.hessian(i, 0, 0))) + t[0] * t[0] * (jit_s[0][1] * (jit_s[1][1] * basis_s.hessian(i, 1, 1) + jit_s[1][0] * basis_s.hessian(i, 1, 0)) + jit_s[0][0] * (jit_s[1][1] * basis_s.hessian(i, 0, 1) + jit_s[1][0] * basis_s.hessian(i, 0, 0)))) * t[1] + (t[1] * t[0] * (jit_s[1][1] * (jit_s[0][1] * basis_s.hessian(i, 1, 1) + jit_s[0][0] * basis_s.hessian(i, 1, 0)) + jit_s[1][0] * (jit_s[0][1] * basis_s.hessian(i, 0, 1) + jit_s[0][0] * basis_s.hessian(i, 0, 0))) + t[0] * t[0] * (jit_s[0][1] * (jit_s[0][1] * basis_s.hessian(i, 1, 1) + jit_s[0][0] * basis_s.hessian(i, 1, 0)) + jit_s[0][0] * (jit_s[0][1] * basis_s.hessian(i, 0, 1) + jit_s[0][0] * basis_s.hessian(i, 0, 0)))) * t[0];

        // Maybe additional corrections of the jump are necessary in below equations by multiplying the flip factor
        double flip_factor = flipped ? 1.0 : -1.0;
        (flipped ? r_s : r_n).accumulate(lfsu_n.child(0), i, E*I*(-sk_dt2un*dtvn_n_0 - dt2vn_n_0 * sk_dtun + beta*sk_dtun*dtvn_n_0));
        (flipped ? r_s : r_n).accumulate(lfsu_n.child(1), i, E*I*(-sk_dt2un*dtvn_n_1 - dt2vn_n_1 * sk_dtun + beta*sk_dtun*dtvn_n_1));

        (flipped ? r_n : r_s).accumulate(lfsu_s.child(0), i, -E*I*(-sk_dt2un*dtvn_s_0 - dt2vn_s_0 * sk_dtun + beta*sk_dtun*dtvn_s_0));
        (flipped ? r_n : r_s).accumulate(lfsu_s.child(1), i, -E*I*(-sk_dt2un*dtvn_s_1 - dt2vn_s_1 * sk_dtun + beta*sk_dtun*dtvn_s_1));
      }
    }
  }

  private:
  double beta;

  const typename GFS::Traits::GridView::IndexSet& is;
  std::map<int, std::vector<std::tuple<std::size_t, double, double>>> element_fibre_intersections;
  std::map<std::pair<int, int>, std::vector<std::tuple<std::size_t, double>>> face_fibre_intersections;
  std::vector<std::shared_ptr<FibreParametrizationBase<2>>> fibre_parametrizations;
  std::vector<double> fibre_modulus;
  std::vector<double> fibre_radii;

  // Store a force vector - code adapted of how the generated code does it.
  // That is not the nicest way of doing it, but we do not need additional
  // template parameters on this class.
  using CoefficientForceLFS = Dune::PDELab::LocalFunctionSpace<GFS>;
  mutable std::shared_ptr<CoefficientForceLFS> coefficient_force_lfs;
  using CoefficientForceVector = Dune::PDELab::Backend::Vector<GFS, double>;
  mutable std::shared_ptr<CoefficientForceVector> coefficient_force_vector;
  mutable std::shared_ptr<const GFS> coefficient_force_gfs;
  using CoefficientForceLFSCache = Dune::PDELab::LFSIndexCache<CoefficientForceLFS>;
  mutable std::shared_ptr<CoefficientForceLFSCache> coefficient_force_lfs_cache;

  using LocalBasisType = typename GFS::template Child<0>::Type::Traits::FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache = Dune::PDELab::LocalBasisCache<LocalBasisType>;
  Cache cache;

  public:
  void setCoefficientForce(std::shared_ptr<const GFS> p_gfs, std::shared_ptr<CoefficientForceVector> p_coefficient_vector)
  {
	coefficient_force_gfs = p_gfs;
	coefficient_force_vector = p_coefficient_vector;
	coefficient_force_lfs = std::make_shared<CoefficientForceLFS>(*coefficient_force_gfs);
	coefficient_force_lfs_cache = std::make_shared<CoefficientForceLFSCache>(*coefficient_force_lfs);
  }
};


template<typename GFS, int dim>
class FibreReinforcedBulkOperator
  : public AbstractLocalOperatorInterface<GFS>
{
  public:
  using BaseOperator = AbstractLocalOperatorInterface<GFS>;
  using BulkOperator = typename OperatorSwitch<GFS, dim>::Elasticity;
  using FibreOperator = EulerBernoulli2DLocalOperator<GFS>;

  using EG = typename BaseOperator::EG;
  using IG = typename BaseOperator::IG;
  using LFSU = typename BaseOperator::LFSU;
  using LFSV = typename BaseOperator::LFSV;
  using X = typename BaseOperator::X;
  using R = typename BaseOperator::R;

  FibreReinforcedBulkOperator(std::shared_ptr<const GFS> gfs,
                              const Dune::ParameterTree& rootparams,
                              const Dune::ParameterTree& params,
                              std::shared_ptr<ElasticMaterialBase<typename GFS::Traits::EntitySet, double>> material)
    : bulkoperator(*gfs, *gfs, params, material)
    , fibreoperator(rootparams, params, gfs)
  {}

  virtual ~FibreReinforcedBulkOperator() = default;

  template<typename CGFS, typename CVEC>
  void setCoefficientForce(std::shared_ptr<CGFS> gfs, std::shared_ptr<CVEC> force)
  {
    bulkoperator.setCoefficientForce(gfs, force);
    fibreoperator.setCoefficientForce(gfs, force);
  }

  void setMaterial(std::shared_ptr<ElasticMaterialBase<typename GFS::Traits::EntitySet, double>> material)
  {
    bulkoperator.setMaterial(material);
  }

  template<typename CGFS, typename CVEC>
  void setCoefficientTraction(std::shared_ptr<CGFS> gfs, std::shared_ptr<CVEC> traction)
  {
    bulkoperator.setCoefficientTraction(gfs, traction);
  }

  virtual void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const override
  {
    bulkoperator.alpha_volume(eg, lfsu, x, lfsv, r);
    fibreoperator.alpha_volume(eg, lfsu, x, lfsv, r);
  }

  virtual void alpha_boundary(const IG& ig, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const override
  {
    bulkoperator.alpha_boundary(ig, lfsu, x, lfsv, r);
  }

  virtual void alpha_skeleton(const IG& ig, const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s, const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, R& r_s, R& r_n) const
  {
    bulkoperator.alpha_skeleton(ig, lfsu_s, x_s, lfsv_s, lfsu_n, x_n, lfsv_n, r_s, r_n);
    fibreoperator.alpha_skeleton(ig, lfsu_s, x_s, lfsv_s, lfsu_n, x_n, lfsv_n, r_s, r_n);
  }

  private:
  typename OperatorSwitch<GFS, dim>::Elasticity bulkoperator;
  EulerBernoulli2DLocalOperator<GFS> fibreoperator;
};


#endif
