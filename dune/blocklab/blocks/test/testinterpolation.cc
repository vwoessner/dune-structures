#include"config.h"

#include<dune/blocklab/blocks/error.hh>
#include<dune/blocklab/blocks/interpolation.hh>
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
  Dune::TestSuite test;
  Dune::ParameterTree config;

  using Grid = Dune::UGGrid<2>;
  using GridProvider = Dune::BlockLab::StructuredSimplexGridProvider<Grid>;
  using VectorProvider = Dune::BlockLab::PkFemVectorProvider<GridProvider, 2>;

  auto grid = std::make_shared<GridProvider>(config);
  auto vector = std::make_shared<VectorProvider>(grid);

  config["solver.blocks"] = "interpolation, error";
  config["interpolation.functions"] = "x * y";
  config["error.analytic"] = "0.0";

  using Context = Dune::BlockLab::ConstructionContext<std::tuple<>, VectorProvider>;
  Context ctx(helper, config, vector);
  ctx.template registerBlock<Dune::BlockLab::InterpolationBlock>("interpolation");
  ctx.template registerBlock<Dune::BlockLab::DiscretizationErrorBlock>("error");

  auto solver = ctx.constructSolver(config.sub("solver"));
  solver->apply();

  // We (ab)use the error node to validate the integral of the interpolate solution
  test.check(std::abs(solver->template param<double>("error") - 1.0/3.0) < 1e-8)
     << "Integration of interpolated function x*y yielded wrong result.";

  return test.exit();
}
