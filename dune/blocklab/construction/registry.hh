#ifndef DUNE_BLOCKLAB_CONSTRUCTION_REGISTRY_HH
#define DUNE_BLOCKLAB_CONSTRUCTION_REGISTRY_HH

/** This header is a curated list of built-in blocks that can
 *  be registered with any construction context.
 */

#include<dune/blocklab/blocks/interpolation.hh>
#include<dune/blocklab/operators/convectiondiffusionfem.hh>

namespace Dune::BlockLab {

  /** Register all the built-in blocks */
  template<typename Context>
  void registerBuiltinBlocks(Context& ctx)
  {
    // Register all the basic blocks from the block subdirectory
    ctx.template registerBlock<InterpolationBlock>("interpolation");

    // Register all the operator-blocks
    ctx.template registerBlock<ConvectionDiffusionFEMBlock>("convectiondiffusionfem");
  }

} // namespace Dune::BlockLab

#endif
