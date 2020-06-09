#include"config.h"

#include<dune/blocklab/grids/concept.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/test/testsuite.hh>

#include<iostream>


class Grid1
{
  public:
  using Grid = int;
  using Parameter = int;
  using EntitySet = int;
};


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;

  test.check(Dune::BlockLab::isGridProvider<Grid1>())
    << "Grid1 failed the concept check but shouldn't";

  return test.exit();
}
