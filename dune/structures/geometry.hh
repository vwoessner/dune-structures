#ifndef DUNE_STRUCTURES_GEOMETRY_HH
#define DUNE_STRUCTURES_GEOMETRY_HH

#include <dune/geometry/multilineargeometry.hh>

#include <memory>
#include <type_traits>
#include <vector>

/// Create a geometry of the displaced configuration
/** This assumes a linear displacement and only operates on the element corners
 */
template<class GridFunction, class Entity, class Geometry>
auto
create_displaced_geometry(
  const std::shared_ptr<GridFunction> displacement_grid_func,
  const Entity& entity,
  const Geometry& geo)
{
  // Check if types are compatible
  using Coordinate = typename Geometry::GlobalCoordinate;
  using RangeType = typename GridFunction::Traits::RangeType;
  static_assert(
    std::is_convertible_v<Coordinate, RangeType>,
    "Displacement field grid function must return compatible coordinates!");

  // Store the displaced corner coordinates
  std::vector<Coordinate> corners_displ;
  for (int i = 0; i < geo.corners(); ++i)
  {
    // Evaluate displacement field at coordinate of original geometry
    RangeType displacement;
    const auto coord = geo.corner(i);
    displacement_grid_func->evaluate(entity, geo.local(coord), displacement);
    corners_displ.push_back(coord + displacement);
  }

  // Create a new geometry and return it
  using DisplacedGeometry = Dune::MultiLinearGeometry<typename Geometry::ctype,
                                                      Geometry::mydimension,
                                                      Geometry::coorddimension>;
  const auto geo_type = geo.type();
  return DisplacedGeometry(geo_type, corners_displ);
}

#endif // DUNE_STRUCTURES_GEOMETRY_HH
