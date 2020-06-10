#include"config.h"

#include<dune/blocklab/blocks/interpolation.hh>
#include<dune/blocklab/construction/context.hh>
#include<dune/blocklab/grids/structured.hh>
#include<dune/blocklab/vectors/pkfem.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertree.hh>
#include<dune/grid/uggrid.hh>

#include<memory>
#include<tuple>

int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  Dune::ParameterTree config;

  using Grid = Dune::UGGrid<2>;
  using GridProvider = Dune::BlockLab::StructuredSimplexGridProvider<Grid>;
  using VectorProvider = Dune::BlockLab::PkFemVectorProvider<GridProvider, 1>;

  auto grid = std::make_shared<GridProvider>(config);
  auto vector = std::make_shared<VectorProvider>(grid);

  config["solver.blocks"] = "interpolation";
  config["interpolation.functions"] = "x + y";

  using Context = Dune::BlockLab::ConstructionContext<std::tuple<>, VectorProvider>;
  Context ctx(helper, config, vector);
  ctx.template registerBlock<Dune::BlockLab::InterpolationBlock>("interpolation");

  auto solver = ctx.constructSolver(config.sub("solver"));
  solver->apply();

  return 0;
}
