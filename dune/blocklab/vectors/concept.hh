#ifndef DUNE_BLOCKLAB_VECTORS_CONCEPT_HH
#define DUNE_BLOCKLAB_VECTORS_CONCEPT_HH

/** A concept for a facility that provides a DoF vector for a BlockSolver
 *
 *  PDELab does not provide a dynamic interface for DoF containers, so
 *  we need to statically provide these. The providing classes need to
 *  fulfill the given concept.
 */

#include<dune/common/concept.hh>

namespace Dune::BlockLab {

  namespace impl {

    struct VectorProviderConcept
    {
      template<typename T>
      auto require(T&& t) -> decltype(
        Dune::Concept::requireType<typename T::FiniteElementMap>()
      );
    };
  } // namespace impl

  template<typename T>
  constexpr bool isVectorProvider()
  {
    return Dune::models<impl::VectorProviderConcept, T>();
  }

} // namespace Dune::BlockLab

#endif
