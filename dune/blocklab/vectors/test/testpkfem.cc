#include"config.h"

#include<dune/blocklab/grids/structured.hh>
#include<dune/blocklab/vectors/pkfem.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/uggrid.hh>

#include<iostream>
#include<memory>


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::ParameterTree config;

  using Grid = Dune::UGGrid<2>;
  using GridProvider = Dune::BlockLab::StructuredSimplexGridProvider<Grid>;
  using VectorProvider = Dune::BlockLab::PkFemVectorProvider<GridProvider, 1>;

  auto grid = std::make_shared<GridProvider>(config);
  auto vector = std::make_shared<VectorProvider>(grid);

  std::cout << vector->getVector()->N() << std::endl;

  return 0;
}
