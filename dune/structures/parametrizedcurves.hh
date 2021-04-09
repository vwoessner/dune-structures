#ifndef DUNE_STRUCTURES_PARAMETERIZEDCURVES_HH
#define DUNE_STRUCTURES_PARAMETERIZEDCURVES_HH

/** A hierarchy of classes to implement parametrized curves in space
 *  as well as an algorithm to intersect those curves with a given grid
 *  and store the relevant data of that intersection.
 */

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <map>
#include <memory>
#include <tuple>

/** An abstract base class for a parametrized curve */
template<int dim>
class FibreParametrizationBase
{
public:
  using GlobalCoordinate = Dune::FieldVector<double, dim>;
  virtual ~FibreParametrizationBase() = default;
  virtual GlobalCoordinate eval(double t) const = 0;
  virtual GlobalCoordinate tangent(double t) const = 0;
};

/** An implementation of a straight line as a parametrized curve
 *
 * This implementation is used in reproduction of the results of the paper
 * 'Cut Finite Element Methods for Linear Elasticity Problems' by Hansbo et al
 * as well as in debugging.
 *
 */
template<int dim>
class StraightFibre : public FibreParametrizationBase<dim>
{
public:
  using GlobalCoordinate = Dune::FieldVector<double, dim>;
  virtual ~StraightFibre() = default;

  StraightFibre(const YAML::Node& param)
    : start(param["start"].as<GlobalCoordinate>())
    , dir(param["end"].as<GlobalCoordinate>() - start)
  {
  }

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

/** A struct that gathers the data of an intersection between a grid and
 *  a parametrized curve
 */
struct FibreGridIntersection
{
  // Intersections with grid cells stored as a mapping from cell index to
  // the parametrization parameter interval of the intersection (encoded as a
  // pair)
  std::map<int, std::pair<double, double>> element_fibre_intersections;

  // Intersections with interior facets of the grid, stored as a mapping of
  // a pair of indices of inside and outside cell to the parametrization
  // parameter of the intersection point.
  std::map<std::pair<int, int>, double> facet_fibre_intersections;
};

/** Find intersections of fibres with the grid. The outline is rougly the
 * following:
 *    - Find the first simplex that the fibre intersects by doing hierarchical
 * search on fibre(0)
 *    - If fibre(0) is outside the grid, do a more expensive method to find the
 * first simplex
 *    - Do a bisection of the fibre interval to find the parameter t where it
 * leaves the simplex
 *    - Store the parameter interval and the simplex index
 *    - Go to the neighboring simplex and repeat until at a boundary
 * intersection or at fibre(1).
 */
template<typename GV>
auto
compute_grid_fibre_intersection(
  GV gv,
  std::shared_ptr<FibreParametrizationBase<GV::dimension>> fibre,
  double tol = 1e-12)
{
  // The return struct that we will fill in the follwoing
  FibreGridIntersection ret;

  const auto& is = gv.indexSet();
  int index = -1;
  auto element = *(gv.template begin<0>());

  // We iterate this algorithm, because in rare scenarios the strategy of going
  // from cell to cell across intersections fails. This is when the curve
  // traverses through a grid vertex from one cell to a non-adjacent one. In
  // these cases we restart the algorithm by doing a grid search for the current
  // cell.
  double t = 0.0;
  while (t + tol < 1.0)
  {
    if (t > tol)
      std::cout << "A restart of the curve cutting algorithm was needed at t="
                << t << std::endl;

    // Find the first cell in the grid that the curve traverses
    for (const auto& e : elements(gv))
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
        for (const auto& e : elements(gv))
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

    // Go from cell to cell by finding the intersection that the curve traverses
    while (index != -1)
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
      ret.element_fibre_intersections[index] = std::make_pair(t, bisect_a);
      t = bisect_a;

      // Find the intersection and the outside cell
      index = -1;
      for (auto intersection : intersections(gv, element))
      {
        if (!intersection.boundary())
        {
          auto geo = intersection.outside().geometry();
          if (referenceElement(geo).checkInside(
                geo.local(fibre->eval(t + tol))))
          {
            int inside_index = is.index(element);
            element = intersection.outside();
            index = is.index(element);
            ret.facet_fibre_intersections[std::make_pair(inside_index, index)] =
              t;
            break;
          }
        }
      }
    }
  }

  return std::move(ret);
}

#endif
