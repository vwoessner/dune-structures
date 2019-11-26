#include"config.h"

#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/uggrid.hh>
#include<dune/pdelab.hh>
#include<dune/structures/material.hh>
#include<dune/structures/onetoone.hh>
#include<dune/structures/transitionsolver.hh>
#include<dune/structures/vonmises.hh>
#include<dune/structures/visualization.hh>
#include<dune/testtools/gridconstruction.hh>

#include<vector>

#include"elasticity_operator.hh"


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // Parse the ini file
  Dune::ParameterTree params;
  Dune::ParameterTreeParser::readINITree("cubeapp.ini", params);

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

  // Set up finite element maps...
  using FEM = Dune::PDELab::PkLocalFiniteElementMap<ES, DF, RangeType, 1>;
  FEM fem(es);

  // Set up grid function spaces...
  using VB = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
  using CASS = Dune::PDELab::ConformingDirichletConstraints;
  using GFS = Dune::PDELab::VectorGridFunctionSpace<ES, FEM, 3, VB, VB, CASS>;
  GFS gfs(es, fem);
  gfs.name("displacement");
  gfs.update();
  std::cout << "Set up a grid function space with " << gfs.size() << " dofs!" << std::endl;

  // Setting up constraints container
  using CC = GFS::ConstraintsContainer<RangeType>::Type;
  CC cc;
  cc.clear();

  // Instantiate the material class
  auto material = parse_material<RangeType>(es, physical, params.sub("material"));

  // Setting up container
  using V = Dune::PDELab::Backend::Vector<GFS, DF>;
  V x(gfs);

  InterpolationTransitionStep<V> interpolation([](auto x) { return 0.0; });
  ConstraintsTransitionStep<V> constraints([](auto x){ return (x[2] < 1e-08) || (x[2] > 1.0 - 1e-8); });
  ElasticitySolverStep<V> elasticity(physical, params);

  double maxdispl = params.get<double>("model.compression");
  ParametrizedTransformationTransitionStep<V, double> trafo(
      [maxdispl](auto u, auto x, double p)
      {
        u[2] = - p * (1.0 - maxdispl) * x[2];
        return u;
      });

  ContinuousVariationTransitionStep<V> compress(10);
  compress.add(trafo);
  compress.add(elasticity);

  VisualizationStep<V> vis;
  SolutionVisualizationStep<V> vissol;
  MPIRankVisualizationStep<V> visrank(helper);
  PhysicalEntityVisualizationStep<V> visphys(physical);
  VonMisesStressVisualizationStep<V> visvm(material);

  vis.add(vissol);
  vis.add(visrank);
  vis.add(visphys);
  vis.add(visvm);

  TransitionSolver<V> solver;
  solver.add(interpolation);
  solver.add(constraints);
  solver.add(compress);
  solver.add(vis);

  solver.apply(x, cc);

  return 0;
}
