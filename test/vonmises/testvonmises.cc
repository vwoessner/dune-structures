#include"config.h"

#include<dune/common/parallel/mpihelper.hh>
#include<dune/pdelab.hh>
#include<dune/structures/material.hh>
#include<dune/structures/vonmises.hh>
#include<dune/testtools/gridconstruction.hh>

#include"linear_elasticity_operator.hh"

int main(int argc, char** argv)
{
	Dune::MPIHelper& mpihelper = Dune::MPIHelper::instance(argc, argv);

  // Parse the ini file
  Dune::ParameterTree params;
  Dune::ParameterTreeParser::readINITree("vonmises.ini", params);

  using RangeType = double;

  // Setup grid (view)...
  using Grid = Dune::UGGrid<3>;
  using GV = Grid::LeafGridView;
  using ES = Dune::PDELab::NonOverlappingEntitySet<GV>;
  using DF = Grid::ctype;
  IniGridFactory<Grid> factory(params.sub("grid"));
  std::shared_ptr<Grid> grid = factory.getGrid();
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
  auto bctype = [&](const auto& is, const auto& xl){ auto x=is.geometry().global(xl); return (abs(x[0]) < 1e-08 ? 1 : 0.0); };
  auto bctype_f = Dune::PDELab::makeBoundaryConditionFromCallable(es, bctype);
  Dune::PDELab::CompositeConstraintsParameters comp_bctype(bctype_f, bctype_f, bctype_f);
  //Dune::PDELab::CompositeConstraintsParameters<decltype(bctype_f), decltype(bctype_f), decltype(bctype_f)> comp_bctype(bctype_f, bctype_f, bctype_f);
  Dune::PDELab::constraints(comp_bctype, gfs, cc);

  // Instantiate the material class
  auto material = std::make_shared<HomogeneousElasticMaterial<ES, double>>(1.25, 1.0);

  // Setting up grid operator
  using LOP = LinearElasticityOperator<GFS, GFS>;
  using MB = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  using GO = Dune::PDELab::GridOperator<GFS, GFS, LOP, MB, DF, RangeType, RangeType, CC, CC>;
  LOP lop(gfs, gfs, params, material);
  MB mb(25);
  GO go(gfs, cc, gfs, cc, lop, mb);

  // Setting up container
  using V = Dune::PDELab::Backend::Vector<GFS, DF>;
  V x(gfs);
  auto dirichlet = [&](const auto& x){ return 0.0; };
  auto dirichlet_f = Dune::PDELab::makeGridFunctionFromCallable(es, dirichlet);
  Dune::PDELab::CompositeGridFunction comp_dirichlet(dirichlet_f, dirichlet_f, dirichlet_f);
  Dune::PDELab::interpolate(comp_dirichlet, gfs, x);

  // Set up the solver...
  using LS = Dune::PDELab::ISTLBackend_SEQ_UMFPack;
  using SLP = Dune::PDELab::StationaryLinearProblemSolver<GO, LS, V>;
  LS ls(false);
  SLP slp(go, ls, x, 1e-12);
  slp.apply();

  // A grid function for the stress
  VonMisesStressGridFunction stress(x, material);

  // Interpolate the stress into a grid function
  using SGFS = Dune::PDELab::GridFunctionSpace<ES, FEM, CASS, VB>;
  SGFS sgfs(es, fem);
  sgfs.name("vonmises");
  using SV = Dune::PDELab::Backend::Vector<SGFS, DF>;
  SV stress_container(sgfs);
  Dune::PDELab::interpolate(stress, sgfs, stress_container);

  // Visualize the stress grid function
  Dune::SubsamplingVTKWriter vtkwriter(gv, Dune::refinementLevels(0));
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, x);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter, sgfs, stress_container);
  vtkwriter.write("vonmises", Dune::VTK::ascii);

  return 0;
}
