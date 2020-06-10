#include"config.h"

#include<dune/blocklab/grids/capabilities.hh>
#include<dune/blocklab/grids/concept.hh>
#include<dune/blocklab/grids/gmsh.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/test/testsuite.hh>
#include<dune/grid/uggrid.hh>
#include<dune/grid/test/gridcheck.hh>


template<typename Provider>
bool checkGrid()
{
  try {
    Dune::ParameterTree p;
    p["filename"] = "unitcube.msh";
    Provider provider(p);
    auto grid = provider.createGrid();
    gridcheck(*grid);
  }
  catch(...)
  {
    return false;
  }
  return true;
}


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;

  using Provider1 = Dune::BlockLab::GMSHGridProvider<Dune::UGGrid<2>>;
  using Provider2 = Dune::BlockLab::GMSHGridProvider<Dune::UGGrid<3>>;

  // Make sure that all providers pass the concept check
  test.check(Dune::BlockLab::isGridProvider<Provider1>())
      << "Provider1 failed the GridProvider concept check";

  test.check(Dune::BlockLab::isGridProvider<Provider2>())
      << "Provider2 failed the GridProvider concept check";

  // Make sure that all capabilities are set correctly
  test.check(Dune::BlockLab::hasSimplices<Provider1>() &&
             Dune::BlockLab::hasCubes<Provider1>() &&
             Dune::BlockLab::isAdaptive<Provider1>())
      << "Provider1 set wrong capabilities";

  test.check(Dune::BlockLab::hasSimplices<Provider2>() &&
             Dune::BlockLab::hasCubes<Provider2>() &&
             Dune::BlockLab::isAdaptive<Provider2>())
      << "Provider2 set wrong capabilities";

  // Make sure that all provider create valid grids
  test.check(checkGrid<Provider1>())
       << "Provider1 created an invalid grid";

  test.check(checkGrid<Provider2>())
       << "Provider2 created an invalid grid";

  return test.exit();
}
