#include"config.h"

#include<dune/blocklab/construction/registry.hh>
#include<dune/blocklab/blocks/test/providersetup.hh>
#include<dune/common/parametertree.hh>
#include<dune/common/parallel/mpihelper.hh>


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  Dune::ParameterTree config;

  // Test a number of vector providers with all blocks added to the context.
  // This is much more of a compile test than anything else

  {
    auto ctx = structured_ug2_p1fem(helper, config);
    registerBuiltinBlocks(ctx);
    ctx.constructSolver(config)->apply();
  }

  {
    auto ctx = structured_ug2_p2fem(helper, config);
    registerBuiltinBlocks(ctx);
    ctx.constructSolver(config)->apply();
  }

  return 0;
}
