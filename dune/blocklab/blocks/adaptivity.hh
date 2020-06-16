#ifndef DUNE_BLOCKLAB_BLOCKS_ADAPTIVITY_HH
#define DUNE_BLOCKLAB_BLOCKS_ADAPTIVITY_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/common/hybridutilities.hh>
#include<dune/common/parametertree.hh>
#include<dune/pdelab/adaptivity/adaptivity.hh>

#include<tuple>

namespace Dune::BlockLab {

  template<typename P, typename V>
  class AdaptivityBlock
    : public ParentBlockBase<P, V>
  {
    public:
    using Traits = BlockTraits<P, V>;

    template<typename Context>
    AdaptivityBlock(Context& ctx, const Dune::ParameterTree& config)
      : ParentBlockBase<P, V>(ctx, config)
    {}

    virtual ~AdaptivityBlock() = default;

    virtual void setup() override
    {
      // Introduce a parameter that notifies steps that the grid has changed.
      // The actual type and value are not particularly relevant: We simply
      // encode that the grid has undergone adaptation. More importantly each
      // update of this parameter can be used as a notification event.
      this->solver->template introduce_parameter<bool>("adapted", false);

      for (auto block : this->blocks)
        block->setup();
    }

    virtual void apply() override
    {
      // TODO: Here, we should check that the grid is unmarked.
      //       All marking should be done by substeps of this one.

      // Recursively call all substeps. These are expected to mark cells for refinement.
      for (auto block : this->blocks)
        block->apply();

      std::cout << "Adapting the grid..." << std::endl;

      /**
       * The following code is horribly complicated. So it might deserve some remarks:
       *
       * * This is necessary because all vectors to be adapted need to be passed to the
       *   `adaptGrid` function in one go.
       * *  PDELab requires a non-const reference to the grid function space, but our way
       *    of not redundantly carrying around gfs's gives us only shared pointers to `const GFS`.
       *    We solve this by casting away the const...
       * *  Another thing that bothers me, where I do not understand PDELab semantics enough:
       *    If two of the vectors in V... use the same grid function space, is it okay to pass
       *    that gfs twice into adapt or do I actually need to group my vectors by gfs object???
       *    I doubt I will ever succeed in doing that in meta-programming, I would much rather
       *    duplicate and invade the PDELab adaptivity interface.
       * *  The `4` below hardcodes the integration order. It should of course be configurable in
       *    some way, but I will not target that before I actually think about coarsening.
       */

      /* This is still broken - I was looking for a way to not expose all vectors through getVectors()
       * but instead iterate over them with an index tuple. Maybe that is just not worth it.
       */
//      auto transfer = std::apply(
//        [](auto... v)
//        {
//          return std::make_tuple(
//            Dune::PDELab::transferSolutions(
//              *const_cast<typename std::remove_const<typename decltype(v)::element_type::GridFunctionSpace>::type*>(v->gridFunctionSpaceStorage().get()),
//              4,
//              *v
//            )...
//          );
//         },
//         this->solver->getVectors()
//       );
//
//      std::apply(
//        [this](auto... v)
//        {
//          Dune::PDELab::adaptGrid(
//            *(this->solver->getGrid()),
//  	  v...
//  	);
//        },
//        transfer
//      );

      // Notify other solver steps that the grid has been adapted
      // They might need to rebuild some data structures.
      this->solver->update_parameter("adapted", true);
    }
  };

} // namespace Dune::BlockLab

#endif
