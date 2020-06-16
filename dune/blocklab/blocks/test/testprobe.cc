#include"config.h"

#include<dune/blocklab/blocks/probe.hh>
#include<dune/blocklab/blocks/interpolation.hh>
#include<dune/blocklab/blocks/test/providersetup.hh>
#include<dune/common/fvector.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertree.hh>
#include<dune/common/test/testsuite.hh>

#include<memory>
#include<tuple>


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;
  Dune::ParameterTree config;

  config["solver.blocks"] = "interpolation, probe";
  config["interpolation.functions"] = "x * y";
  config["probe.position"] = "0.5 0.5";
  config["probe.name"] = "myprobe";

  auto ctx = structured_ug2_p2fem(helper, config);
  ctx.template registerBlock<Dune::BlockLab::InterpolationBlock>("interpolation");
  ctx.template registerBlock<Dune::BlockLab::ProbeBlock>("probe");

  auto solver = ctx.constructSolver(config.sub("solver"));
  solver->apply();

  // We (ab)use the error node to validate the integral of the interpolate solution
  test.check(std::abs(solver->template param<Dune::FieldVector<double, 1>>("myprobe") - 0.25) < 1e-8)
     << "Probing function x*y at (0.5, 0.5) yielded wrong result.";

  return test.exit();
}
