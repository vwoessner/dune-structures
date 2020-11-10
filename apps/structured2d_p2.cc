#include"config.h"

/** This is the main app that allows to access all of BlockLab's core
 *  capabilities through one executable. Of course, it takes horribly
 *  long to compile. For downstream projects, it makes much more sense
 *  to provide their own app, that follows the implementation patterns
 *  used in this app, but only includes the necessary grid provider,
 *  vector providers and blocks.
 */

#include<dune/blocklab/app.hh>
#include<dune/blocklab/events.hh>
#include<dune/blocklab/init.hh>
#include<dune/blocklab/grids.hh>
#include<dune/blocklab/vectors.hh>

#include<dune/blocklab/utilities/tuple.hh>
#include<dune/pdelab/common/partitionviewentityset.hh>

#include<dune/structures/gridprovider.hh>
#include<dune/structures/material.hh>
#include<dune/structures/registry.hh>

#include<memory>
#include<tuple>


int main(int argc, char** argv)
{
  auto init = Dune::BlockLab::initBlockLab(argc, argv);

  using GridProvider = Dune::BlockLab::StructuredSimplexGridProvider<Dune::UGGrid<2>>;
  auto grid = [](const auto& c)
    { 
      return std::make_shared<GridProvider>(c);
    };

  auto vectors = std::make_tuple(
    std::make_tuple("Displacement Field", "P2 Lagrange Element", [](auto gp)
      {
        using GridProvider = typename decltype(gp)::element_type;
        auto leaf = std::make_shared<Dune::BlockLab::PkFemVectorProvider<GridProvider, 2>>(gp);
        return Dune::BlockLab::fieldProvider<GridProvider::Grid::dimension>(leaf);
      }),
    std::make_tuple("Body Force", "P2 Lagrange Element", [](auto gp)
      {
        using GridProvider = typename decltype(gp)::element_type;
        auto leaf = std::make_shared<Dune::BlockLab::PkFemVectorProvider<GridProvider, 2>>(gp);
        return Dune::BlockLab::fieldProvider<GridProvider::Grid::dimension>(leaf);
      }),
    std::make_tuple("Traction Force", "P2 Lagrange Element", [](auto gp)
      {
        using GridProvider = typename decltype(gp)::element_type;
        auto leaf = std::make_shared<Dune::BlockLab::PkFemVectorProvider<GridProvider, 2>>(gp);
        return Dune::BlockLab::fieldProvider<GridProvider::Grid::dimension>(leaf);
      })
  );

  using Material = std::shared_ptr<ElasticMaterialBase<Dune::PDELab::OverlappingEntitySet<typename GridProvider::Grid::Traits::LeafGridView>, double>>;
  using ParameterTuple = std::tuple<std::shared_ptr<std::vector<int>>, Material>;

  Dune::BlockLab::BlockLabApp app(init, grid, vectors, registerStructuresBlocks, ParameterTuple{});

  app.addDefaultRunner();
  app.addFrontendExporter();
  app.addHelpMessage();
  app.setTitle("Structured 2D Cube Grid - Quadratic Elements");

  app.run();

  return 0;
}
