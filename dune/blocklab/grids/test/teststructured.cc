#include"config.h"

#include<dune/blocklab/grids/concept.hh>
#include<dune/blocklab/grids/structured.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/test/testsuite.hh>
#include<dune/grid/uggrid.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/test/gridcheck.hh>


template<typename Provider>
bool checkGrid()
{
  try {
    Dune::ParameterTree p;

    if (Provider::dim == 1)
      p["N"] = "1";
    if (Provider::dim == 2)
      p["N"] = "1 1";
    if (Provider::dim == 3)
      p["N"] = "1 1 1";

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

  using Provider1 = Dune::BlockLab::StructuredSimplexGridProvider<Dune::UGGrid<2>>;
  using Provider2 = Dune::BlockLab::StructuredCubeGridProvider<Dune::UGGrid<2>>;
  using Provider3 = Dune::BlockLab::StructuredSimplexGridProvider<Dune::UGGrid<3>>;
  using Provider4 = Dune::BlockLab::StructuredCubeGridProvider<Dune::UGGrid<3>>;
  using Provider5 = Dune::BlockLab::StructuredCubeGridProvider<Dune::YaspGrid<2>>;
  using Provider6 = Dune::BlockLab::StructuredCubeGridProvider<Dune::YaspGrid<3>>;

  // Make sure that all providers pass the concept check
  test.check(Dune::BlockLab::isGridProvider<Provider1>())
      << "Provider1 failed the GridProvider concept check";

  test.check(Dune::BlockLab::isGridProvider<Provider2>())
      << "Provider2 failed the GridProvider concept check";

  test.check(Dune::BlockLab::isGridProvider<Provider3>())
      << "Provider3 failed the GridProvider concept check";

  test.check(Dune::BlockLab::isGridProvider<Provider4>())
      << "Provider4 failed the GridProvider concept check";

  test.check(Dune::BlockLab::isGridProvider<Provider5>())
      << "Provider5 failed the GridProvider concept check";

  test.check(Dune::BlockLab::isGridProvider<Provider6>())
      << "Provider6 failed the GridProvider concept check";

  test.check(checkGrid<Provider1>())
      << "Provider1 created an invalid grid";

  test.check(checkGrid<Provider2>())
      << "Provider2 created an invalid grid";

  test.check(checkGrid<Provider3>())
      << "Provider3 created an invalid grid";

  test.check(checkGrid<Provider4>())
      << "Provider4 created an invalid grid";

  test.check(checkGrid<Provider5>())
      << "Provider5 created an invalid grid";

  test.check(checkGrid<Provider6>())
      << "Provider6 created an invalid grid";

  return test.exit();
}
