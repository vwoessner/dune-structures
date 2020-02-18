#ifndef DUNE_STRUCTURES_GRIDCONSTRUCTION_HH
#define DUNE_STRUCTURES_GRIDCONSTRUCTION_HH

#include<dune/common/filledarray.hh>
#include<dune/common/parametertree.hh>
#include<dune/grid/uggrid.hh>
#include<dune/pdelab.hh>
#include<dune/structures/gmshfactory.hh>

#include<cstdlib>
#include<fstream>
#include<memory>
#include<tuple>


template<int dim=3>
auto construct_grid(Dune::MPIHelper& helper, const Dune::ParameterTree& params, char** argv)
{
  using GridType = Dune::UGGrid<dim>;
  using GV = typename GridType::LeafGridView;
  //using ES = Dune::PDELab::NonOverlappingEntitySet<GV>;
  using ES = Dune::PDELab::OverlappingEntitySet<GV>;
  using DF = typename GridType::ctype;

  auto type = params.get<std::string>("type", "structured");
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

  if (type != "structured")
    DUNE_THROW(Dune::Exception, "Your grid type is not known!");

  auto ll = params.get<Dune::FieldVector<double, dim>>("lowerleft", Dune::FieldVector<double, dim>(0.0));
  auto ur = params.get<Dune::FieldVector<double, dim>>("upperright", Dune::FieldVector<double, dim>(1.0));
  auto N = params.get<std::array<unsigned int, dim>>("N", Dune::filledArray<dim, unsigned int>(10));

  auto grid = Dune::StructuredGridFactory<GridType>::createSimplexGrid(ll, ur, N);
  grid->loadBalance();
  auto physical = std::make_shared<std::vector<int>>(grid->size(0), 0);
  GV gv = grid->leafGridView();
  ES es(grid->leafGridView());

  return std::make_tuple(std::shared_ptr<GridType>(std::move(grid)), es, physical);
}

#endif
