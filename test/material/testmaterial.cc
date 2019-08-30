#include"config.h"

#include<dune/common/parametertree.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/uggrid.hh>
#include<dune/pdelab.hh>
#include<dune/structures/gmshfactory.hh>
#include<dune/structures/material.hh>

#include<iostream>

// The generated operator
#include"material_test_operator.hh"

int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // Parse the ini file
  Dune::ParameterTree config;
  Dune::ParameterTreeParser::readINITree("testmaterial.ini", config);

  // Construct the grid
  using GridType = Dune::UGGrid<3>;
  PhysicalEntityGmshFactory<GridType> factory(helper, config.sub("grid"));
  auto grid = factory.getGrid();
  auto physical = factory.getPhysical();

  // Parse material information
  using GV = GridType::LeafGridView;
  using ES = Dune::PDELab::NonOverlappingEntitySet<GV>;
  ES es(grid->leafGridView());
  auto material = parse_material<double>(es, physical, config.sub("material"));

  // Build a local operator
  using FEM = Dune::PDELab::PkLocalFiniteElementMap<ES, double, double, 1>;
  using CASS = Dune::PDELab::NoConstraints;
  using VB = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
  using GFS = Dune::PDELab::GridFunctionSpace<ES, FEM, CASS, VB>;
  using LOP = MaterialTestOperator<GFS, GFS>;
  FEM fem(es);
  GFS gfs(es, fem);
  LOP lop(gfs, gfs, config, material);

  // Build a grid operator from that
  using MB = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  using CC = GFS::ConstraintsContainer<double>::Type;
  using GO = Dune::PDELab::GridOperator<GFS, GFS, LOP, MB, double, double, double, CC, CC>;
  MB mb(25);
  CC cc;
  GO go(gfs, cc, gfs, cc, lop, mb);

  // Calculate a residual
  using V = Dune::PDELab::Backend::Vector<GFS, double>;
  V x(gfs, 0.0), r(gfs, 0.0);
  go.residual(x, r);
}
