#ifndef DUNE_BLOCKLAB_BLOCKS_TEST_PROVIDERSETUP_HH
#define DUNE_BLOCKLAB_BLOCKS_TEST_PROVIDERSETUP_HH

/** This header defines some test scenarios to be reused in
 *  unit tests for blocks to reduce the amount of code duplication.
 */

#include<dune/blocklab/construction/context.hh>
#include<dune/blocklab/grids/structured.hh>
#include<dune/blocklab/vectors/pkfem.hh>
#include<dune/common/parametertree.hh>
#include<dune/grid/uggrid.hh>

#include<memory>
#include<tuple>


auto structured_ug2_p1fem(Dune::MPIHelper& helper, const Dune::ParameterTree& config)
{
  using Grid = Dune::UGGrid<2>;
  using GridProvider = Dune::BlockLab::StructuredSimplexGridProvider<Grid>;
  using VectorProvider = Dune::BlockLab::PkFemVectorProvider<GridProvider, 1>;

  auto grid = std::make_shared<GridProvider>(config);
  auto vector = std::make_shared<VectorProvider>(grid);

  using Context = Dune::BlockLab::ConstructionContext<std::tuple<>, VectorProvider>;
  Context ctx(helper, config, vector);

  return ctx;
}


auto structured_ug2_p2fem(Dune::MPIHelper& helper, const Dune::ParameterTree& config)
{
  using Grid = Dune::UGGrid<2>;
  using GridProvider = Dune::BlockLab::StructuredSimplexGridProvider<Grid>;
  using VectorProvider = Dune::BlockLab::PkFemVectorProvider<GridProvider, 2>;

  auto grid = std::make_shared<GridProvider>(config);
  auto vector = std::make_shared<VectorProvider>(grid);

  using Context = Dune::BlockLab::ConstructionContext<std::tuple<>, VectorProvider>;
  Context ctx(helper, config, vector);

  return ctx;
}

#endif
