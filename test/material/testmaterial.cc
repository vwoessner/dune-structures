#include"config.h"

#include<dune/common/parametertree.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/uggrid.hh>
#include<dune/structures/material.hh>

#include<iostream>


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // Parse the ini file
  Dune::ParameterTree config;
  Dune::ParameterTreeParser::readINITree("material.ini", config);

  // Construct the grid on rank 0 to later distribute it
  using GridType = Dune::UGGrid<3>;
  Dune::GridFactory<GridType> factory;

  // Physical entity information container
  std::vector<int> boundary, entity;

  if (helper.rank() == 0)
  {
    Dune::GmshReader<GridType>::read(factory, "testmaterial.msh", boundary, entity, true, false);
  }
  auto grid = std::shared_ptr<GridType>(factory.createGrid());
  grid->loadBalance();

  using GV = GridType::LeafGridView;
  const GV& gv = grid->leafGridView();

  auto material = parse_material<double>(gv, entity, config);

  for (auto e : Dune::elements(gv))
    std::cout << material.first_lame(e, Dune::FieldVector<double, 3>(0.0)) << " "
              << material.second_lame(e, Dune::FieldVector<double, 3>(0.0)) << std::endl;
}
