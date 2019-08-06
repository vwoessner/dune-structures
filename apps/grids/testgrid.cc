#include"config.h"

#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/test/gridcheck.hh>
#include<dune/grid/uggrid.hh>

#include<memory>


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // Construct the grid on rank 0 to later distribute it
  using GridType = Dune::UGGrid<3>;
  Dune::GridFactory<GridType> factory;
  if (helper.rank() == 0)
  {
    Dune::GmshReader<GridType>::read(factory, "testgrid.msh", true, false);
  }
  auto grid = std::shared_ptr<GridType>(factory.createGrid());

  // This is only feasible for small meshes!
  // gridcheck(*grid);

  return 0;
}
