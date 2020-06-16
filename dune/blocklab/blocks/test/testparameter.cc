#include"config.h"

#include<dune/blocklab/blocks/parameter.hh>
#include<dune/blocklab/blocks/test/providersetup.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertree.hh>
#include<dune/common/test/testsuite.hh>


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;
  Dune::ParameterTree config;

  config["solver.blocks"] = "parameter";
  config["parameter.datatype"] = "int";
  config["parameter.value"] = "42";
  config["parameter.name"] = "myparam";

  auto ctx = structured_ug2_p1fem(helper, config);
  ctx.template registerBlock<Dune::BlockLab::ParameterBlock>("parameter");

  auto solver = ctx.constructSolver(config.sub("solver"));
  solver->apply();

  test.check(solver->template param<int>("myparam") == 42)
     << "Retrieval of int parameter myparam failed";

  return test.exit();
}
