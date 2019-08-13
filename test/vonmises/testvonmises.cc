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
  using Grid = Dune::YaspGrid<3, Dune::EquidistantCoordinates<RangeType, 3>>;
  using GV = Grid::LeafGridView;
  using DF = Grid::ctype;
  IniGridFactory<Grid> factory(params);
  std::shared_ptr<Grid> grid = factory.getGrid();
  GV gv = grid->leafGridView();

  // Set up finite element maps...
  using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV, DF, RangeType, 1>;
  FEM fem(gv);

  // Set up grid function spaces...
  using VB = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
  using CASS = Dune::PDELab::ConformingDirichletConstraints;
  using GFS = Dune::PDELab::VectorGridFunctionSpace<GV, FEM, 3, VB, VB, CASS>;
  GFS gfs(gv, fem);
  gfs.name("displacement");

  // Setting up container
  using V = Dune::PDELab::Backend::Vector<GFS, DF>;
  V x(gfs);

  // A grid function for the stress
  auto material = std::make_shared<HomogeneousElasticMaterial<GV, double>>(1.0, 1.0);
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
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter, sgfs, stress_container);
  vtkwriter.write("vonmises", Dune::VTK::ascii);

  return 0;
}
