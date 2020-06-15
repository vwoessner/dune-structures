#include"config.h"

#include<dune/blocklab/blocks/constraints.hh>
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
  Dune::TestSuite test;
  Dune::ParameterTree config;

  config["solver.blocks"] = "constraints";
  config["constraints.functions"] = "x < 1e-8";

  auto ctx = structured_ug2_p1fem(helper, config);
  ctx.template registerBlock<Dune::BlockLab::ConstraintsBlock>("constraints");

  auto solver = ctx.constructSolver(config.sub("solver"));
  solver->apply();

  // We (ab)use the error node to validate the integral of the interpolate solution
  test.check(solver->template getConstraintsContainer<0>()->size() == 11)
     << "Number of constrained DoFs for x<eps yielded wrong result.";

  return test.exit();
}
