#ifndef DUNE_BLOCKLAB_GRIDS_CAPABILITIES_HH
#define DUNE_BLOCKLAB_GRIDS_CAPABILITIES_HH

/** This header provides some capabilities for grid providers.
 *  The has* and is* methods can be used like concepts in order
 *  to restrict functionality to matching grid providers at compile
 *  time. The Has* and Is* structs from the namespace Capabilities
 *  need to be specialized by all grid providers. Note that, while
 *  similar, these capabilities differ from those in dune-grid:
 *  E.g. a provider may set hasCubes to false for Dune::UGGrid although
 *  UGGrid supports cubes if the provider uses UGGrid as a simplex grid.
 */

#include<dune/grid/uggrid.hh>


namespace Dune::BlockLab {

  namespace Capabilities {

    template<typename T>
    struct HasSimplices
    {
      static constexpr bool value = false;
    };

    template<typename T>
    struct HasCubes
    {
      static constexpr bool value = false;
    };

    template<typename T>
    struct IsAdaptive
    {
      static constexpr bool value = false;
    };

    // Some specializations based on the grid type
    template<template<typename> typename Provider, int dim>
    struct IsAdaptive<Provider<Dune::UGGrid<dim>>>
    {
      static constexpr bool value = true;
    };

  }

  template<typename T>
  constexpr bool hasSimplices()
  {
    return Capabilities::HasSimplices<T>::value;
  }

  template<typename T>
  constexpr bool hasCubes()
  {
    return Capabilities::HasCubes<T>::value;
  }

  template<typename T>
  constexpr bool isAdaptive()
  {
    return Capabilities::IsAdaptive<T>::value;
  }

  template<typename T>
  constexpr bool hasPrismsPyramids()
  {
    return hasCubes<T>() && isAdaptive<T>();
  }

  template<typename T>
  constexpr bool hasOneGeometry()
  {
    if constexpr (hasCubes<T>() && isAdaptive<T>())
      return false;
    if constexpr (hasCubes<T>() && hasSimplices<T>())
      return false;

    return true;
  }

}

#endif
