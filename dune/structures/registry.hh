#ifndef DUNE_STRUCTURES_REGISTRY_HH
#define DUNE_STRUCTURES_REGISTRY_HH

#include <dune/blocklab/construction/registry.hh>
#include <dune/structures/elasticity.hh>
#include <dune/structures/material.hh>
#include <dune/structures/onetoone.hh>
#include <dune/structures/reinforcedoperator.hh>
#include <dune/structures/visualization.hh>

auto registerStructuresBlocks = [](auto& ctx) {
  Dune::BlockLab::registerBuiltinBlocks(ctx);

  // Add the structures-specific blocks to the solver
  ctx.template registerBlock<ElasticityOperatorBlock>("elasticity");
  ctx.template registerBlock<MaterialInitializationBlock>("material");
  ctx.template registerBlock<FibreReinforcedElasticityOperatorBlock>(
    "reinforced_operator");
  ctx.template registerBlock<ElasticityMassOperatorBlock>("mass_operator");
  ctx.template registerBlock<OneToOneMappingCheckerBlock>("onetoone");
  ctx.template registerBlock<FibreDistanceVisualizationBlock>(
    "vis_fibredistance");
  ctx.template registerBlock<PhysicalEntityVisualizationBlock>("vis_physical");
  ctx.template registerBlock<VonMisesStressVisualizationBlock>("vis_vonmises");
  ctx.template registerBlock<StressEVVisualizationBlock>("vis_stress");
};

#endif
