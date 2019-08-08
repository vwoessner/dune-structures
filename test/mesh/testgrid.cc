#include"config.h"

#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/test/gridcheck.hh>
#include<dune/grid/uggrid.hh>

#include<map>
#include<memory>


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
    Dune::GmshReader<GridType>::read(factory, "testgrid.msh", boundary, entity, true, false);
  }
  auto grid = std::shared_ptr<GridType>(factory.createGrid());

  // Print a summary of encountered physical entity groups
  std::map<int, int> occurences;
  for(auto v: entity)
	  ++occurences[v];

  std::cout << "\nEncountered physical entity groups:" << std::endl;
  for(auto v: occurences)
	 std::cout << "Group " << v.first << ": " << v.second << " times" << std::endl;

  // This is only feasible for small meshes!
  // gridcheck(*grid);

  return 0;
}
