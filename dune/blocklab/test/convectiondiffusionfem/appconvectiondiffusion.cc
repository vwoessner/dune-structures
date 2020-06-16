#include"config.h"

#include<dune/blocklab/construction/context.hh>
#include<dune/blocklab/construction/registry.hh>
#include<dune/blocklab/grids/structured.hh>
#include<dune/blocklab/vectors/pkfem.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertree.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/test/testsuite.hh>
#include<dune/grid/uggrid.hh>

#include<memory>


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  std::string inifile;
  if (argc == 1)
    inifile = "convectiondiffusionfem.ini";
  else
    inifile = argv[1];

  // Parse the ini file
  Dune::ParameterTree config;
  Dune::ParameterTreeParser::readOptions(argc, argv, config);
  Dune::ParameterTreeParser::readINITree(inifile, config, false);

  using Grid = Dune::UGGrid<2>;
  using GridProvider = Dune::BlockLab::StructuredSimplexGridProvider<Grid>;
  using VectorProvider = Dune::BlockLab::PkFemVectorProvider<GridProvider, 1>;

  auto grid = std::make_shared<GridProvider>(config.sub("grid"));
  auto vector = std::make_shared<VectorProvider>(grid);

  using Context = Dune::BlockLab::ConstructionContext<std::tuple<>, VectorProvider>;
  Context ctx(helper, config, vector);
  Dune::BlockLab::registerBuiltinBlocks(ctx);

  auto solver = ctx.constructSolver(config.sub("solver"));
  solver->apply();

  Dune::TestSuite test;
  // The error threshold here is the L2-error of the interpolant of the solution
  test.check(solver->template param<double>("l2error") < 3e-5);

  return test.exit();
}
