#include"config.h"

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertree.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/structures/elasticity.hh>
#include<dune/structures/solverconstruction.hh>


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  using RangeType = double;

  // Parse the ini file
  Dune::ParameterTree params;
  Dune::ParameterTreeParser::readINITree("solver1.ini", params);

  // Construct the grid and loadbalance it
  using GridType = Dune::UGGrid<3>;
  using GV = GridType::LeafGridView;
  using ES = Dune::PDELab::NonOverlappingEntitySet<GV>;
  using DF = GridType::ctype;

  Dune::FieldVector<double, 3> ll(0.0);
  Dune::FieldVector<double, 3> ur(1.0);
  auto N = params.get<std::array<unsigned int, 3>>("grid.N", {1, 1, 1});

  auto grid = Dune::StructuredGridFactory<GridType>::createSimplexGrid(ll, ur, N);
  auto physical = std::make_shared<std::vector<int>>(grid->size(0), 0);
  GV gv = grid->leafGridView();
  ES es(grid->leafGridView());

  auto [x, cc] = elasticity_setup(es);
  using V = std::remove_reference<decltype(*x)>::type;

  ConstructionContext<V> ctx(helper, params, es, physical);
  auto solver = ctx.construct(params);

  solver->apply(x, cc);

  return 0;
}
