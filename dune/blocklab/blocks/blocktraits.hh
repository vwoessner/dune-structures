#ifndef DUNE_BLOCKLAB_BLOCKS_BLOCKTRAITS_HH
#define DUNE_BLOCKLAB_BLOCKS_BLOCKTRAITS_HH

#include<dune/common/fvector.hh>
#include<dune/common/parametertree.hh>

#include<string>
#include<tuple>


namespace Dune::BlockLab {

  // The forward declaration of the BlockSolver class
  template<typename P, typename V>
  class BlockSolver;

  template<typename P, typename V, std::size_t i=0>
  class BlockTraits;

  template<std::size_t i, typename... P, typename... V>
  class BlockTraits<std::tuple<P...>, std::tuple<V...>, i>
  {
    public:
    // The transition solver for this
    using Solver = BlockSolver<std::tuple<P...>, std::tuple<V...>>;

    // Export the vector and grid dependent types
    using Vector = typename std::tuple_element<i, std::tuple<V...>>::type;
    using GridFunctionSpace = typename Vector::GridFunctionSpace;
    using GridView = typename GridFunctionSpace::Traits::GridViewType;
    static constexpr int dim = GridView::dimension;
    using Grid = typename GridView::Traits::Grid;
    using EntitySet = typename GridFunctionSpace::Traits::EntitySet;
    using ctype = typename Grid::ctype;
    using Entity = typename GridView::template Codim<0>::Entity;
    using GlobalCoordinate = typename Entity::Geometry::GlobalCoordinate;
    using Range = typename Vector::field_type;
    using ConstraintsContainer = typename GridFunctionSpace::template ConstraintsContainer<Range>::Type;
    using VectorBackend = typename GridFunctionSpace::Traits::Backend;
    using Parameter = typename Solver::Parameter;
  };

} // namespace Dune::BlockLab

#endif
