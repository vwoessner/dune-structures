#include"config.h"

#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/uggrid.hh>
#include<dune/pdelab.hh>
#include<dune/structures/material.hh>
#include<dune/structures/vonmises.hh>
#include<dune/testtools/gridconstruction.hh>

#include<vector>

#include"linear_elasticity_operator.hh"

int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // Parse the ini file
  Dune::ParameterTree params;
  Dune::ParameterTreeParser::readINITree("simplecell.ini", params);

  using RangeType = double;

  // Construct the grid on rank 0 to later distribute it
  using GridType = Dune::UGGrid<3>;
  using GV = GridType::LeafGridView;
  using DF = GridType::ctype;
  Dune::GridFactory<GridType> factory;

  // Physical entity information container
  std::vector<int> boundary, entity;

  if (helper.rank() == 0)
  {
    Dune::GmshReader<GridType>::read(factory, "simplecell.msh", boundary, entity, true, false);
  }
  auto grid = std::shared_ptr<GridType>(factory.createGrid());
  GV gv = grid->leafGridView();

  // Set up finite element maps...
  using FEM = Dune::PDELab::PkLocalFiniteElementMap<GV, DF, RangeType, 1>;
  FEM fem(gv);

  // Set up grid function spaces...
  using VB = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
  using CASS = Dune::PDELab::ConformingDirichletConstraints;
  using GFS = Dune::PDELab::VectorGridFunctionSpace<GV, FEM, 3, VB, VB, CASS>;
  GFS gfs(gv, fem);
  gfs.name("displacement");
  gfs.update();
  std::cout << "Set up a grid function space with " << gfs.size() << " dofs!" << std::endl;

  // Setting up constraints container
  using CC = GFS::ConstraintsContainer<RangeType>::Type;
  CC cc;
  cc.clear();
  auto bctype = [&](const auto& x){ return (x[0] < -0.5 + 1e-08 ? 1 : 0); };
  auto bctype_f = Dune::PDELab::makeBoundaryConditionFromCallable(gv, bctype);
  Dune::PDELab::CompositeConstraintsParameters comp_bctype(bctype_f, bctype_f, bctype_f);
  Dune::PDELab::constraints(comp_bctype, gfs, cc);
  std::cout << "Set up a constraints container with " << cc.size() << " dofs!" << std::endl;

  // Instantiate the material class
  auto material = parse_material<RangeType>(gv, entity, params.sub("material"));
//  auto material = std::make_shared<HomogeneousElasticMaterial<GV, double>>(1.25, 1.0);

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
  auto dirichlet_f = Dune::PDELab::makeGridFunctionFromCallable(gv, dirichlet);
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
  using SGFS = Dune::PDELab::GridFunctionSpace<GV, FEM, CASS, VB>;
  SGFS sgfs(gv, fem);
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
