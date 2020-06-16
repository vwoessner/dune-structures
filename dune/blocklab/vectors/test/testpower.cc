#include"config.h"

#include<dune/blocklab/grids/structured.hh>
#include<dune/blocklab/vectors/concept.hh>
#include<dune/blocklab/vectors/pkfem.hh>
#include<dune/blocklab/vectors/power.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/test/testsuite.hh>
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
  auto pvector = Dune::BlockLab::powerProvider<5>("power", std::make_shared<VectorProvider>(grid));

  Dune::TestSuite test;

  // Check the concept definition
  test.check(Dune::BlockLab::isVectorProvider<Dune::BlockLab::PowerVectorProvider<VectorProvider, 5>>())
    << "PowerVectorProvider failed the concept check!";

  // Assert that the number of DoFs is indeed the correct multiple compared to leaf space
  test.check(vector->getVector()->N() * 5 == pvector->getVector()->N())
    << "Numer of DoFs in power space is not multiple.";

  return test.exit();
}
