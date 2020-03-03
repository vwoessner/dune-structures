#include"config.h"

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertree.hh>
#include<dune/structures/elasticity.hh>
#include<dune/structures/gridconstruction.hh>
#include<dune/structures/solverconstruction.hh>

#include<memory>


template<int dim, int degree>
void apply(Dune::MPIHelper& helper, const Dune::ParameterTree& params, char** argv)
{
  auto [grid, es, physical] = construct_grid<dim>(helper, params.sub("grid"), argv);
  auto [x, cc] = elasticity_setup<degree>(es);
  using V = typename std::remove_reference<decltype(*x)>::type;

  ConstructionContext<V> ctx(helper, params, es, physical);
  ctx.setVectors(x);
  ctx.setConstraintsContainers(cc);


  ctx.template registerOperator<0>("elasticity",
      [](const auto& ctx, const auto& p)
      {
        auto vec = ctx.template getVector<0>();
        auto gfs = vec->gridFunctionSpaceStorage();
        return std::make_shared<typename OperatorSwitch<typename V::GridFunctionSpace, dim>::Elasticity>
            (*gfs, *gfs, p, nullptr);
      });

  auto solver = ctx.construct(params.sub("solver"));

  solver->apply();
}


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // Parse the ini file
  Dune::ParameterTree params;
  Dune::ParameterTreeParser::readINITree(argv[1], params);

  // Dispatch on grid dimension
  auto dim = params.get<int>("grid.dimension", 3);
  auto degree = params.get<int>("solver.degree", 1);

  if (dim == 2)
  {
    if (degree == 1)
      apply<2, 1>(helper, params, argv);
    if (degree == 2)
      apply<2, 2>(helper, params, argv);
  }
  if (dim == 3)
  {
    if (degree == 1)
      apply<3, 1>(helper, params, argv);
    if (degree == 2)
      apply<3, 2>(helper, params, argv);
  }
  return 0;
}
