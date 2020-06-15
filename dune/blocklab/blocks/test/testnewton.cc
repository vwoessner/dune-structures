#include"config.h"

#include<dune/blocklab/blocks/newton.hh>
#include<dune/blocklab/blocks/test/providersetup.hh>
#include<dune/blocklab/operators/convectiondiffusionfem.hh>
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

  config["solver.blocks"] = "convectiondiffusionfem, newton";
  config["newton.operator"] = "convectiondiffusionfem_operator";

  auto ctx = structured_ug2_p1fem(helper, config);
  ctx.template registerBlock<Dune::BlockLab::NewtonSolverBlock>("newton");
  ctx.template registerBlock<Dune::BlockLab::ConvectionDiffusionFEMBlock>("convectiondiffusionfem");

  auto solver = ctx.constructSolver(config.sub("solver"));
  solver->apply();

  return test.exit();
}
