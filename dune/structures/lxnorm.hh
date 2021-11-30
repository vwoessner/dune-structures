#ifndef DUNE_STRUCTURES_LXNORM_HH
#define DUNE_STRUCTURES_LXNORM_HH

#include <cmath>

#include <dune/grid/common/gridview.hh>
#include <dune/pdelab/common/quadraturerules.hh>

#include <dune/structures/geometry.hh>

/// Calculate the norm of a function over a displaced geometry
/**
 *  \param u Function to compute the norm on
 *  \param displacement_grid_func Grid function representing the displacement
 *    field
 *  \param order Order of the norm
 *  \param quadrature_order Order of the quadrature rule for integrating
 */
template<typename U, class GridFunction>
typename U::Traits::RangeFieldType
lxnorm_trf(const U& u,
           const std::shared_ptr<GridFunction> displacement_grid_func,
           const int order = 2,
           const int quadrature_order = 1)
{
  using namespace std;
  using RF = typename U::Traits::RangeFieldType;
  RF sum = 0.0;
  for (auto& elem : Dune::elements(u.getGridView()))
  {
    // Use a transformed geometry
    const auto geo = elem.geometry();
    const auto geo_trf =
      create_displaced_geometry(displacement_grid_func, elem, geo);
    for (const auto& it : Dune::quadratureRule(geo_trf, quadrature_order))
    {
      // evaluate the given grid functions at integration point
      typename U::Traits::RangeType u_val;
      u.evaluate(elem, it.position(), u_val);

      // accumulate error
      const auto val = std::accumulate(
        begin(u_val), end(u_val), RF(0.0), [order](const auto a, const auto b) {
          return a + pow(b, order);
        });
      sum += val * it.weight() * geo_trf.integrationElement(it.position());
    }
  }
  return pow(sum, 1.0 / order);
}

// Calculate the L2 norm of a function
// template<typename U>
// auto
// lxnorm(const U& u, const int order = 2, const int quadrature_order = 1)
// {
//   using namespace std;
//   return pow(lxnormx(u, order, quadrature_order), 1 / order);
// }

#endif // DUNE_STRUCTURES_LXNORM_HH
