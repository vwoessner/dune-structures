#include"config.h"

#include<dune/blocklab/blocks/error.hh>
#include<dune/blocklab/construction/context.hh>
#include<dune/blocklab/grids/structured.hh>
#include<dune/blocklab/vectors/pkfem.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertree.hh>
#include<dune/common/test/testsuite.hh>
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

  config["solver.blocks"] = "error";
  config["error.analytic"] = "1.0";

  using Context = Dune::BlockLab::ConstructionContext<std::tuple<>, VectorProvider>;
  Context ctx(helper, config, vector);
  ctx.template registerBlock<Dune::BlockLab::DiscretizationErrorBlock>("error");

  auto solver = ctx.constructSolver(config.sub("solver"));
  solver->apply();

  Dune::TestSuite test;
  test.check(std::abs(solver->template param<double>("error") - 1.0) < 1e-8)
    << "Integration of constant 1 in error context failed";

  return test.exit();
}
