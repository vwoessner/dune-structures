#ifndef DUNE_BLOCKLAB_GRIDS_CONCEPT_HH
#define DUNE_BLOCKLAB_GRIDS_CONCEPT_HH

/** A concept for a facility that provides a grid for a BlockSolver
 *
 *  Due to the different grid implementations having different grid types,
 *  a dynamic interface for these factories is not possible. We solve the
 *  issue by defining a concept instead!
 */

#include<dune/common/concept.hh>

namespace Dune::BlockLab {

  namespace impl {

    struct GridProviderConcept
    {
      template<class T>
      auto require(T&& t) -> decltype(
	Dune::Concept::requireType<typename T::Grid>(),
	Dune::Concept::requireType<typename T::Parameter>(),
	Dune::Concept::requireType<typename T::EntitySet>()
      );
    };

  } // namespace impl

  template<typename T>
  constexpr bool isGridProvider()
  {
    return Dune::models<impl::GridProviderConcept, T>();
  }

} // namespace Dune::BlockLab

#endif
