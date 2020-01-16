#ifndef DUNE_STRUCTURES_GRIDCONSTRUCTION_HH
#define DUNE_STRUCTURES_GRIDCONSTRUCTION_HH

#include<dune/common/parametertree.hh>
#include<dune/grid/uggrid.hh>
#include<dune/pdelab.hh>
#include<dune/structures/gmshfactory.hh>

#include<cstdlib>
#include<fstream>
#include<memory>
#include<tuple>


auto construct_grid(Dune::MPIHelper& helper, const Dune::ParameterTree& params, char** argv)
{
  using GridType = Dune::UGGrid<3>;
  using GV = GridType::LeafGridView;
  using ES = Dune::PDELab::NonOverlappingEntitySet<GV>;
  using DF = GridType::ctype;

  auto type = params.get<std::string>("type", "structured");
  if (type == "structured")
  {
    auto ll = params.get<Dune::FieldVector<double, 3>>("lowerleft", Dune::FieldVector<double, 3>(0.0));
    auto ur = params.get<Dune::FieldVector<double, 3>>("lowerleft", Dune::FieldVector<double, 3>(1.0));
    auto N = params.get<std::array<unsigned int, 3>>("N", {10, 10, 10});

    auto grid = Dune::StructuredGridFactory<GridType>::createSimplexGrid(ll, ur, N);
    grid->loadBalance();
    auto physical = std::make_shared<std::vector<int>>(grid->size(0), 0);
    GV gv = grid->leafGridView();
    ES es(grid->leafGridView());

    return std::make_tuple(std::shared_ptr<GridType>(std::move(grid)), es, physical);
  }
  if (type == "cell")
  {
    // Trigger python mesh generator
    using namespace std::literals;
    auto command = std::string(VIRTUALENV_WRAPPER_DIR) + "/run-in-dune-env generate_cell_mesh "s + std::string(argv[1]) + " >meshgen.log";
    std::system(command.c_str());

    // Read the generated mesh file
    PhysicalEntityGmshFactory<GridType> factory(helper, params);
    auto grid = factory.getGrid();
    auto physical = factory.getPhysical();
    GV gv = grid->leafGridView();
    ES es(grid->leafGridView());

    return std::make_tuple(grid, es, physical);
  }
}

#endif
