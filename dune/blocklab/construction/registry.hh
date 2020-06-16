#ifndef DUNE_BLOCKLAB_CONSTRUCTION_REGISTRY_HH
#define DUNE_BLOCKLAB_CONSTRUCTION_REGISTRY_HH

/** This header is a curated list of built-in blocks that can
 *  be registered with any construction context.
 */

#include<dune/blocklab/blocks/constraints.hh>
#include<dune/blocklab/blocks/control.hh>
#include<dune/blocklab/blocks/error.hh>
#include<dune/blocklab/blocks/filelogger.hh>
#include<dune/blocklab/blocks/instationary.hh>
#include<dune/blocklab/blocks/interpolation.hh>
#include<dune/blocklab/blocks/linearsolver.hh>
#include<dune/blocklab/blocks/newton.hh>
#include<dune/blocklab/blocks/parameter.hh>
#include<dune/blocklab/blocks/probe.hh>
#include<dune/blocklab/blocks/variation.hh>
#include<dune/blocklab/blocks/visualization.hh>
#include<dune/blocklab/operators/convectiondiffusionfem.hh>


namespace Dune::BlockLab {

  /** Register all the built-in blocks */
  template<typename Context>
  void registerBuiltinBlocks(Context& ctx)
  {
    // Register all the basic blocks from the block subdirectory
    ctx.template registerBlock<ConstraintsBlock>("constraints");
    ctx.template registerBlock<ContinuousVariationBlock>("continuousvariation");
    ctx.template registerBlock<DiscreteVariationBlock>("discretevariation");
    ctx.template registerBlock<DiscretizationErrorBlock>("error");
    ctx.template registerBlock<FileLoggerBlock>("filelogger");
    ctx.template registerBlock<InterpolationBlock>("interpolation");
    ctx.template registerBlock<LinearSolverBlock>("linearsolver");
    ctx.template registerBlock<NewtonSolverBlock>("newton");
    ctx.template registerBlock<ProbeBlock>("probe");
    ctx.template registerBlock<ParameterBlock>("parameter");
    ctx.template registerBlock<RepeatBlock>("repeat");
    ctx.template registerBlock<TimestepperBlock>("timestepper");
    ctx.template registerBlock<VisualizationBlock>("visualization");
    ctx.template registerBlock<IndexSetVisualizationBlock>("vis_indexset");
    ctx.template registerBlock<MPIRankVisualizationBlock>("vis_mpirank");
    ctx.template registerBlock<VectorVisualizationBlock>("vis_vector");

    // Register all the operator-blocks
    ctx.template registerBlock<ConvectionDiffusionFEMBlock>("convectiondiffusionfem");
  }

} // namespace Dune::BlockLab

#endif
