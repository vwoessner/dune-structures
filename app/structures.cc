#include"config.h"

/** This is the main app that allows to access all of BlockLab's core
 *  capabilities through one executable. Of course, it takes horribly
 *  long to compile. For downstream projects, it makes much more sense
 *  to provide their own app, that follows the implementation patterns
 *  used in this app, but only includes the necessary grid provider,
 *  vector providers and blocks.
 */

#include<dune/blocklab/construction/combinatorics.hh>
#include<dune/blocklab/grids.hh>
#include<dune/blocklab/utilities/tuple.hh>
#include<dune/blocklab/vectors.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertree.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/pdelab/common/partitionviewentityset.hh>

#include<dune/structures/elasticity.hh>
#include<dune/structures/gridprovider.hh>
#include<dune/structures/material.hh>
#include<dune/structures/onetoone.hh>
#include<dune/structures/visualization.hh>
#include<dune/structures/reinforcedoperator.hh>

#include<memory>
#include<tuple>


template<typename... P>
struct ParameterGatherer
{
  using Materials = std::tuple<std::shared_ptr<ElasticMaterialBase<Dune::PDELab::OverlappingEntitySet<typename P::Grid::Traits::LeafGridView>, double>>...>;
};


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // Parse the ini file
  Dune::ParameterTree config;
  Dune::ParameterTreeParser::readOptions(argc, argv, config);
  Dune::ParameterTreeParser::readINITree(argv[1], config, false);

  auto grid_providers = std::make_tuple(
//    std::make_pair("cell3d", [](const auto& c){ return std::make_shared<StructuresGridProvider<3>>(c); })//,
    std::make_pair("structuredug2d", [](const auto& c){ return std::make_shared<Dune::BlockLab::StructuredSimplexGridProvider<Dune::UGGrid<2>>>(c); })//,
//    std::make_pair("structuredug3d", [](const auto& c){ return std::make_shared<Dune::BlockLab::StructuredSimplexGridProvider<Dune::UGGrid<3>>>(c); })
  );

  auto vector_providers = std::make_tuple(
//    std::make_pair("p1fem",
//                   [](const auto& c, auto gp){
//                     using GridProvider = typename decltype(gp)::element_type;
//                     auto leaf = std::make_shared<Dune::BlockLab::PkFemVectorProvider<GridProvider, 1>>(gp);
//                     return Dune::BlockLab::fieldProvider<GridProvider::Grid::dimension>(leaf);
//                   })//,
    std::make_pair("p2fem",
                   [](const auto& c, auto gp){
                     using GridProvider = typename decltype(gp)::element_type;
                     auto leaf = std::make_shared<Dune::BlockLab::PkFemVectorProvider<GridProvider, 2>>(gp);
                     return Dune::BlockLab::fieldProvider<GridProvider::Grid::dimension>(leaf);
                   })
  );

  // The registration function that we are using
  auto reg = [](auto& ctx){
    Dune::BlockLab::registerBuiltinBlocks(ctx);

    // Add the structures-specific blocks to the solver
    ctx.template registerBlock<ElasticityOperatorBlock>("elasticity_operator");
    ctx.template registerBlock<MaterialInitializationBlock>("material");
    ctx.template registerBlock<FibreReinforcedElasticityOperatorBlock>("reinforced_operator");
    ctx.template registerBlock<ElasticityMassOperatorBlock>("mass_operator");
    ctx.template registerBlock<OneToOneMappingCheckerBlock>("onetoone");
    ctx.template registerBlock<FibreDistanceVisualizationBlock>("vis_fibredistance");
    ctx.template registerBlock<PhysicalEntityVisualizationBlock>("vis_physical");
    ctx.template registerBlock<VonMisesStressVisualizationBlock>("vis_vonmises");
  };

  // Construct additional types to inject into the parameter system
  using MaterialTuple = ParameterGatherer<
      StructuresGridProvider<3>,
      Dune::BlockLab::StructuredSimplexGridProvider<Dune::UGGrid<2>>,
      Dune::BlockLab::StructuredSimplexGridProvider<Dune::UGGrid<3>>
      >::Materials;

  using ParameterTuple = Dune::BlockLab::tuple_cat_t<std::tuple<std::shared_ptr<std::vector<int>>>, MaterialTuple>;

  Dune::BlockLab::instantiate_combinatorics<3, ParameterTuple>(helper, config, grid_providers, vector_providers, reg);

  return 0;
}
