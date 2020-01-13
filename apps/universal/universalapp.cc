#include"config.h"

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertree.hh>
#include<dune/structures/elasticity.hh>
#include<dune/structures/gridconstruction.hh>
#include<dune/structures/solverconstruction.hh>


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // Parse the ini file
  Dune::ParameterTree params;
  Dune::ParameterTreeParser::readINITree(argv[1], params);

  auto [grid, es, physical] = construct_grid(helper, params.sub("grid"), argv);
  auto [x, cc] = elasticity_setup(es);
  using V = std::remove_reference<decltype(*x)>::type;

  ConstructionContext<V> ctx(helper, params, es, physical);
  auto solver = ctx.construct(params.sub("solver"));

  solver->apply(x, cc);

  return 0;
}
