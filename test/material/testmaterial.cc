#include"config.h"

#include<dune/common/parametertree.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/uggrid.hh>
#include<dune/pdelab.hh>
#include<dune/structures/material.hh>

#include<iostream>

// The generated operator
#include"material_test_operator.hh"

int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // Parse the ini file
  Dune::ParameterTree config;
  Dune::ParameterTreeParser::readINITree("material.ini", config);

  // Construct the grid
  using GridType = Dune::UGGrid<3>;
  Dune::GridFactory<GridType> factory;
  std::vector<int> boundary, entity;
  if (helper.rank() == 0)
  {
    Dune::GmshReader<GridType>::read(factory, "testmaterial.msh", boundary, entity, true, false);
  }
  auto grid = std::shared_ptr<GridType>(factory.createGrid());
  grid->loadBalance();

  // Parse material information
  using GV = GridType::LeafGridView;
  const GV& gv = grid->leafGridView();
  auto material = parse_material<double>(gv, entity, config);

  // Build a local operator
  using FEM = Dune::PDELab::PkLocalFiniteElementMap<GV, double, double, 1>;
  using CASS = Dune::PDELab::NoConstraints;
  using VB = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
  using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM, CASS, VB>;
  using LOP = MaterialTestOperator<GFS, GFS, decltype(material)>;

  FEM fem(gv);
  GFS gfs(gv, fem);
  LOP lop(gfs, gfs, config, material);
}
