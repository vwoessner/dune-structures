#include"config.h"

#include<dune/blocklab/grids/structured.hh>
#include<dune/blocklab/vectors/concept.hh>
#include<dune/blocklab/vectors/pkfem.hh>
#include<dune/blocklab/vectors/composite.hh>
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
  using VectorProvider1 = Dune::BlockLab::PkFemVectorProvider<GridProvider, 1>;
  using VectorProvider2 = Dune::BlockLab::PkFemVectorProvider<GridProvider, 2>;

  auto grid = std::make_shared<GridProvider>(config);
  auto vector1 = std::make_shared<VectorProvider1>(grid);
  auto vector2 = std::make_shared<VectorProvider2>(grid);
  auto cvector = Dune::BlockLab::compositeProvider("both",
						   std::make_shared<VectorProvider1>(grid),
						   std::make_shared<VectorProvider2>(grid));

  Dune::TestSuite test;

  // Check the concept definition
  test.check(Dune::BlockLab::isVectorProvider<Dune::BlockLab::CompositeVectorProvider<VectorProvider1, VectorProvider2>>())
    << "PowerVectorProvider failed the concept check!";

  // Assert that the number of DoFs is indeed the sum of DoFs in the leaf spaces
  test.check(vector1->getVector()->N() + vector2->getVector()->N() == cvector->getVector()->N())
    << "Numer of DoFs in composite space is not sum of children.";

  return test.exit();
}
