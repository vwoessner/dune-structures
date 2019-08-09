#include"config.h"

#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/uggrid.hh>
#include<dune/structures/material.hh>

#include<iostream>


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

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

  HomogeneousElasticMaterial<GV, double> mat1(1.0, 1.0);

  MaterialCollection<GV, double> coll(gv, entity);
  coll.add_material(0, mat1);
  coll.add_material(1, std::make_shared<HomogeneousElasticMaterial<GV, double>>(2.0, 2.0));

  std::cout << coll.first_lame(*gv.begin<0>(), Dune::FieldVector<double, 3>(0.0)) << std::endl;
  std::cout << coll.second_lame(*gv.begin<0>(), Dune::FieldVector<double, 3>(0.0)) << std::endl;
}
