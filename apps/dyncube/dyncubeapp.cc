#include"config.h"

#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/uggrid.hh>
#include<dune/pdelab.hh>
#include<dune/structures/elasticity.hh>
#include<dune/structures/material.hh>
#include<dune/structures/onetoone.hh>
#include<dune/structures/timecapsule.hh>
#include<dune/structures/transitionsolver.hh>
#include<dune/structures/vonmises.hh>
#include<dune/structures/visualization.hh>
#include<dune/testtools/gridconstruction.hh>

#include<vector>


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // Parse the ini file
  Dune::ParameterTree params;
  Dune::ParameterTreeParser::readINITree(argv[1], params);

  using RangeType = double;

  // Construct the grid and loadbalance it
  using GridType = Dune::UGGrid<3>;
  using GV = GridType::LeafGridView;
  using ES = Dune::PDELab::NonOverlappingEntitySet<GV>;
  using DF = GridType::ctype;

  Dune::FieldVector<double, 3> ll(0.0);
  Dune::FieldVector<double, 3> ur(1.0);
  auto N = params.get<std::array<unsigned int, 3>>("grid.N", {10, 10, 10});
  std::cout << "N: " << N[0] << " " << N[1] << " " << N[2] << std::endl;

  auto grid = Dune::StructuredGridFactory<GridType>::createSimplexGrid(ll, ur, N);
  auto physical = std::make_shared<std::vector<int>>(grid->size(0), 0);
  GV gv = grid->leafGridView();
  ES es(grid->leafGridView());

  auto [x, cc] = elastodynamics_setup(es);
  using V = std::remove_reference<decltype(*x)>::type;

  TransitionSolver<V> solver;
  TimeCapsule tc;

  // Instantiate the material class
  auto material = parse_material<RangeType>(es, physical, params.sub("material"));

  // Interpolation of initial condition
  auto zero = [](auto x) { return 0.0; };
  InterpolationTransitionStep<V> interpolation(zero);
  solver.add(interpolation);

  // Set up constraints container
  auto clamp = [](auto x){ return (x[2] < 1e-08) || (x[2] > 1.0 - 1e-8); };
  ConstraintsTransitionStep<V> constraints(clamp);
  solver.add(constraints);

  // Set up visualization
  VisualizationStep<V, true> vis;
  SolutionVisualizationStep<V, true> vissol;
  vis.add(vissol);

  // The time stepper
  InstationarySolverStep<V> instat(0.0, 1.0, 0.1);
  solver.add(instat);

  // Shifting the bottom to the right
  ElastoDynamicsSolverStep<V> elastodyn(es, physical, params, tc,
                                        [&tc](auto x){ return 0.2 * x[2] * tc.getTime(); }, zero, zero,
                                        [](auto x){ return x[2] * 0.2;}, zero, zero);
  instat.add(elastodyn);

  // Set up released constraints container
  auto rclamp = [](auto x){ return x[2] < 1e-08; };
  ConstraintsTransitionStep<V> release_constraints(rclamp);
  solver.add(release_constraints);

  // The time stepper
  InstationarySolverStep<V> instat2(1.0, 1.3, 0.001);
  solver.add(instat2);

  // Shifting the bottom to the right
  ElastoDynamicsSolverStep<V> elastodyn2(es, physical, params);
  instat2.add(elastodyn2);
  instat2.add(vis);

  // Run the thing
  solver.apply(x, cc);

  return 0;
}
