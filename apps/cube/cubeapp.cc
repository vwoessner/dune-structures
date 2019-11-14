#include"config.h"

#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/uggrid.hh>
#include<dune/pdelab.hh>
#include<dune/structures/material.hh>
#include<dune/structures/onetoone.hh>
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
  Dune::ParameterTreeParser::readINITree("nonlinear.ini", params);

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

  // Setting up grid operator
  using LOP = ElasticityOperator<GFS, GFS>;
  LOP lop(gfs, gfs, params, material);

  // Setting up container
  using V = Dune::PDELab::Backend::Vector<GFS, DF>;
  V x(gfs);

  auto compression = params.get<double>("model.compression");
  InterpolationTransitionStep<V> interpolation([&](const auto& x){ return (compression - 1.0) * x[2]; },
                                               [](auto x) { return 0.0; },
                                               [](auto x) { return 0.0; });
  ConstraintsSolverStep<V> constraints([](auto x){ return (x[2] < 1e-08) || (x[2] > 1.0 - 1e-8); });
  NewtonSolverStep<V, LOP> newton(lop);

  auto interpolation = std::make_shared<InterpolationTransitionStep<V>>([&](const auto& x){ return (compression - 1.0) * x[2]; },
                                                 [](auto x) { return 0.0; },
                                                 [](auto x) { return 0.0; });
  auto constraints = std::make_shared<ConstraintsSolverStep<V>>([](auto x){ return (x[2] < 1e-08) || (x[2] > 1.0 - 1e-8); });
  auto newton = std::make_shared<NewtonSolverStep<V, LOP>>(lop);

  TransitionSolver<V> solver;
  solver.add(interpolation);
  solver.add(constraints);
  solver.add(newton);

  solver.apply(x, cc);

  // A grid function for the stress
  VonMisesStressGridFunction stress(x, material);

  // Interpolate the stress into a grid function
  using P0FEM = Dune::PDELab::P0LocalFiniteElementMap<DF, RangeType, 3>;
  P0FEM p0fem(Dune::GeometryTypes::simplex(3));
  using P0GFS = Dune::PDELab::GridFunctionSpace<ES, P0FEM, CASS, VB>;
  P0GFS p0gfs(es, p0fem);
  p0gfs.name("vonmises");
  using SV = Dune::PDELab::Backend::Vector<P0GFS, DF>;
  SV stress_container(p0gfs);
  Dune::PDELab::interpolate(stress, p0gfs, stress_container);

  // Visualize the stress grid function
  Dune::VTKWriter vtkwriter(es.gridView());
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, x);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter, p0gfs, stress_container);
  vtkwriter.addCellData(*physical, "gmshPhysical");
  write_rankdata(vtkwriter, helper, es.gridView());
  vtkwriter.write("output", Dune::VTK::ascii);

  if (!is_onetoone(gf))
  {
    is_onetoone(gf, true);
    std::cout << "By the way, your solution is garbage!" << std::endl;
  }
  return 0;
}
