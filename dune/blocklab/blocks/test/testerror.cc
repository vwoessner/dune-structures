#include"config.h"

#include<dune/blocklab/blocks/error.hh>
#include<dune/blocklab/blocks/test/providersetup.hh>
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

  config["solver.blocks"] = "error";
  config["error.analytic"] = "1.0";

  auto ctx = structured_ug2_p1fem(helper, config);
  ctx.template registerBlock<Dune::BlockLab::DiscretizationErrorBlock>("error");

  auto solver = ctx.constructSolver(config.sub("solver"));
  solver->apply();

  Dune::TestSuite test;
  test.check(std::abs(solver->template param<double>("error") - 1.0) < 1e-8)
    << "Integration of constant 1 in error context failed";

  return test.exit();
}
