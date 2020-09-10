#ifndef DUNE_STRUCTURES_EULERBERNOULLI_HH
#define DUNE_STRUCTURES_EULERBERNOULLI_HH

/** Local operators to implement the beam-truss elements for
 *  Euler-Bernoulli beams used in
 *
 *  Hansbo P., Larson M. and Larsson K.
 *  Cut Finite Element Methods for Linear Elasticity Problems
 */

#include<dune/blocklab/operators/virtualinterface.hh>
#include<dune/blocklab/utilities/enumerate.hh>
#include<dune/blocklab/utilities/stringsplit.hh>
#include<dune/geometry/affinegeometry.hh>
#include<dune/pdelab.hh>
#include<dune/structures/elasticity.hh>
#include<dune/structures/material.hh>
#include<dune/structures/parametrizedcurves.hh>

#include<algorithm>
#include<array>
#include<map>
#include<tuple>


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
template<typename GFS, typename FGFS>
class EulerBernoulli2DLocalOperator
  :  public Dune::PDELab::FullSkeletonPattern
  ,  public Dune::PDELab::FullVolumePattern
  ,  public Dune::PDELab::LocalOperatorDefaultFlags
  ,  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  ,  public Dune::PDELab::NumericalJacobianVolume<EulerBernoulli2DLocalOperator<GFS, FGFS>>
  ,  public Dune::PDELab::NumericalJacobianSkeleton<EulerBernoulli2DLocalOperator<GFS, FGFS>>
{
  public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  enum { doPatternSkeleton = true };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doAlphaSkeleton  = true };

  EulerBernoulli2DLocalOperator(const YAML::Node& rootparams, const YAML::Node& params, std::shared_ptr<const GFS> gfs)
    : Dune::PDELab::NumericalJacobianVolume<EulerBernoulli2DLocalOperator<GFS, FGFS>>(1e-9)
    , Dune::PDELab::NumericalJacobianSkeleton<EulerBernoulli2DLocalOperator<GFS, FGFS>>(1e-9)
    , gv(gfs->gridView())
  {
    // Some debugging switches to reduce recompilations in debugging
    bool verbose = params["verbose"].as<bool>(false);
    enable_volume = params["enable_volume"].as<bool>(true);
    enable_skeleton = params["enable_skeleton"].as<bool>(true);

	// Extract stabilization parameter
    beta = params["stabilization_parameter"].as<double>(1000.0);

    // Parse fibres from the configuration
    for (auto fibre: rootparams["grid"]["fibres"])
    {
      if (fibre["shape"].as<std::string>() == "cylinder")
        fibre_parametrizations.push_back(std::make_shared<StraightFibre<2>>(fibre));
      else
        DUNE_THROW(Dune::Exception, "Fibre shape not supported!");

      // Parse the Youngs modulus of the fibre.
      auto modulus = fibre["youngs_modulus"].as<double>(0.0);
      if (modulus == 0.0)
        std::cout << "youngs_modulus not set for all fibres!" << std::endl;
      fibre_modulus.push_back(modulus);

      // Parse the fibre radii
      auto radius = fibre["radius"].as<double>(0.0);
      if (radius == 0.0)
        std::cout << "radius not set for all fibres!" << std::endl;
      fibre_radii.push_back(radius);
    }
    std::cout << "Parsed a total of " << fibre_parametrizations.size() << " fibres from the configuration." << std::endl;

    compute_grid_intersection();

    if (verbose)
    {
      std::cout << "Fibre intersection summary:" << std::endl;
      for (auto [fibindex, intersection] : Dune::BlockLab::enumerate(fibre_intersections))
      {
        for (auto [cellindex, range] : intersection.element_fibre_intersections)
        {
          auto [start, end] = range;
          std::cout << "Cell " << cellindex << " intersects fibre " << fibindex
                    << " on the curve interval [" << start << "," << end << "]" << std::endl;
        }

        for (auto [indexpair, tparam] : intersection.facet_fibre_intersections)
        {
          auto [inside, outside] = indexpair;
          std::cout << "Facet between cell " << inside << " and " << outside << " intersects fibre "
                    << fibindex << " at t=" << tparam << std::endl;
        }
      }
    }
  }

  void compute_grid_intersection()
  {
    fibre_intersections.resize(fibre_parametrizations.size());
    std::transform(fibre_parametrizations.begin(),
		   fibre_parametrizations.end(),
		   fibre_intersections.begin(),
		   [this](auto f){ return compute_grid_fibre_intersection(this->gv, f); });
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    if (!enable_volume)
      return;

    // See whether something needs to be done on this cell
    auto entity = eg.entity();
    const auto& is = gv.indexSet();

    for (std::size_t fibindex=0; fibindex<fibre_parametrizations.size(); ++fibindex)
    {
      // Check whether this fiber actually intersects the cell
      auto fibintersection = fibre_intersections[fibindex];
      auto it = fibintersection.element_fibre_intersections.find(is.index(entity));
      if (it == fibintersection.element_fibre_intersections.end())
        continue;

      auto fibre = fibre_parametrizations[fibindex];

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

      auto start = fibre->eval(it->second.first);
      auto stop = fibre->eval(it->second.second);

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
        auto t = fibre->tangent(it->second.first + ip.position() * (it->second.second - it->second.first));

        // Evaluate the body force vector
        Dune::FieldVector<double, 2> force(0.0);
        for (std::size_t k = 0; k < 2; ++k)
          for (std::size_t i = 0; i < child_0.size(); ++i)
            force[k] += local_coefficient_force_vector(lfsv.child(k), i) * basis.function(i);

        // Extract physical parameter of the fibre
        auto E = fibre_modulus[fibindex];
        auto d = fibre_radii[fibindex];
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
    if (!enable_skeleton)
      return;

    // See whether something needs to be done on this intersection
    const auto& is = gv.indexSet();

    for (std::size_t fibindex=0; fibindex<fibre_parametrizations.size(); ++fibindex)
    {
      // The notion of inside and outside is quite complex in this situation.
      // We are traversing the 2D grid, whose intersections have inside and outside.
      // But these do not necessarily match the notion of inside and outside of
      // the 1D discretization. We therefore need to remap these references in some
      // cases.
      bool flipped = false;

      // Check whether this fiber actually intersects the cell
      auto fibintersection = fibre_intersections[fibindex];
      auto it = fibintersection.facet_fibre_intersections.find(std::make_pair(is.index(ig.inside()), is.index(ig.outside())));
      if (it == fibintersection.facet_fibre_intersections.end())
      {
        it = fibintersection.facet_fibre_intersections.find(std::make_pair(is.index(ig.outside()), is.index(ig.inside())));
        if (it == fibintersection.facet_fibre_intersections.end())
          continue;

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

      auto fibre = fibre_parametrizations[fibindex];
      auto tparam = fibre->eval(it->second);

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
      auto t = fibre->tangent(it->second);

      // Extract physical parameter of the fibre
      auto E = fibre_modulus[fibindex];
      auto d = fibre_radii[fibindex];
      auto A = d;
      auto I = (d*d*d) / 12.0;

      // Compute the penalty factor
      auto penalty = beta / ig.geometry().volume();

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

        double flip_factor = flipped ? 1.0 : -1.0;
        (flipped ? r_s : r_n).accumulate(lfsu_n.child(0), i, E*I*(-sk_dt2un*dtvn_n_0 - flip_factor * dt2vn_n_0 * sk_dtun + flip_factor * penalty * sk_dtun * dtvn_n_0));
        (flipped ? r_s : r_n).accumulate(lfsu_n.child(1), i, E*I*(-sk_dt2un*dtvn_n_1 - flip_factor * dt2vn_n_1 * sk_dtun + flip_factor * penalty * sk_dtun * dtvn_n_1));

        (flipped ? r_n : r_s).accumulate(lfsu_s.child(0), i, -E*I*(-sk_dt2un*dtvn_s_0 - flip_factor * dt2vn_s_0 * sk_dtun + flip_factor * penalty * sk_dtun * dtvn_s_0));
        (flipped ? r_n : r_s).accumulate(lfsu_s.child(1), i, -E*I*(-sk_dt2un*dtvn_s_1 - flip_factor * dt2vn_s_1 * sk_dtun + flip_factor * penalty * sk_dtun * dtvn_s_1));
      }
    }
  }

  private:

  bool enable_volume;
  bool enable_skeleton;

  double beta;

  typename GFS::Traits::GridView gv;
  std::vector<std::shared_ptr<FibreParametrizationBase<2>>> fibre_parametrizations;
  std::vector<FibreGridIntersection> fibre_intersections;
  std::vector<double> fibre_modulus;
  std::vector<double> fibre_radii;

  // Store a force vector - code adapted of how the generated code does it.
  // That is not the nicest way of doing it, but we do not need additional
  // template parameters on this class.
  using CoefficientForceLFS = Dune::PDELab::LocalFunctionSpace<FGFS>;
  mutable std::shared_ptr<CoefficientForceLFS> coefficient_force_lfs;
  using CoefficientForceVector = Dune::PDELab::Backend::Vector<FGFS, double>;
  mutable std::shared_ptr<CoefficientForceVector> coefficient_force_vector;
  mutable std::shared_ptr<const FGFS> coefficient_force_gfs;
  using CoefficientForceLFSCache = Dune::PDELab::LFSIndexCache<CoefficientForceLFS>;
  mutable std::shared_ptr<CoefficientForceLFSCache> coefficient_force_lfs_cache;

  using LocalBasisType = typename GFS::template Child<0>::Type::Traits::FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache = Dune::PDELab::LocalBasisCache<LocalBasisType>;
  Cache cache;

  public:
  void setCoefficientForce(std::shared_ptr<const FGFS> p_gfs, std::shared_ptr<CoefficientForceVector> p_coefficient_vector)
  {
	coefficient_force_gfs = p_gfs;
	coefficient_force_vector = p_coefficient_vector;
	coefficient_force_lfs = std::make_shared<CoefficientForceLFS>(*coefficient_force_gfs);
	coefficient_force_lfs_cache = std::make_shared<CoefficientForceLFSCache>(*coefficient_force_lfs);
  }
};


template<typename GFS, typename FGFS, typename TGFS, int dim>
class FibreReinforcedBulkOperator
  : public Dune::BlockLab::AbstractLocalOperatorInterface<GFS>
{
  public:
  using BaseOperator = Dune::BlockLab::AbstractLocalOperatorInterface<GFS>;
  using BulkOperator = typename OperatorSwitch<GFS, dim, FGFS, TGFS>::Elasticity;
  using FibreOperator = EulerBernoulli2DLocalOperator<GFS, FGFS>;

  using EG = typename BaseOperator::EG;
  using IG = typename BaseOperator::IG;
  using LFSU = typename BaseOperator::LFSU;
  using LFSV = typename BaseOperator::LFSV;
  using X = typename BaseOperator::X;
  using R = typename BaseOperator::R;

  FibreReinforcedBulkOperator(std::shared_ptr<const GFS> gfs,
                              const YAML::Node& rootparams,
                              const YAML::Node& params,
                              std::shared_ptr<ElasticMaterialBase<typename GFS::Traits::EntitySet, double>> material)
    : bulkoperator(*gfs, *gfs, Dune::ParameterTree(), material)
    , fibreoperator(rootparams, params, gfs)
  {}

  virtual ~FibreReinforcedBulkOperator() = default;

  void compute_grid_intersection()
  {
    fibreoperator.compute_grid_intersection();
  }

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
  typename OperatorSwitch<GFS, dim, FGFS, TGFS>::Elasticity bulkoperator;
  EulerBernoulli2DLocalOperator<GFS, FGFS> fibreoperator;
};


#endif
