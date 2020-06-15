#include"config.h"

#include<dune/blocklab/blocks/error.hh>
#include<dune/blocklab/blocks/interpolation.hh>
#include<dune/blocklab/blocks/test/providersetup.hh>
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

  config["solver.blocks"] = "interpolation, error";
  config["interpolation.functions"] = "x * y";
  config["error.analytic"] = "0.0";

  auto ctx = structured_ug2_p2fem(helper, config);
  ctx.template registerBlock<Dune::BlockLab::InterpolationBlock>("interpolation");
  ctx.template registerBlock<Dune::BlockLab::DiscretizationErrorBlock>("error");

  auto solver = ctx.constructSolver(config.sub("solver"));
  solver->apply();

  // We (ab)use the error node to validate the integral of the interpolate solution
  test.check(std::abs(solver->template param<double>("error") - 1.0/3.0) < 1e-8)
     << "Integration of interpolated function x*y yielded wrong result.";

  return test.exit();
}
