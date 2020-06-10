#include"config.h"

/** The unit test for the ConstructionContext object.
 * It tests the registration of blocks and their construction
 * with dummy blocks.
 */

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/construction/context.hh>
#include<dune/blocklab/grids/structured.hh>
#include<dune/blocklab/vectors/pkfem.hh>
#include<dune/common/parametertree.hh>
#include<dune/grid/uggrid.hh>
#include<dune/pdelab.hh>

#include<iostream>
#include<tuple>


// A standard block that does not specify the vector index
template<typename P, typename V>
class DummyBlock1
  : public Dune::BlockLab::BlockBase<P, V>
{
  public:
  template<typename Context>
  DummyBlock1(Context&, const Dune::ParameterTree&)
  {
    std::cout << "Block Dummy1 created!" << std::endl;
  }

  virtual ~DummyBlock1() = default;

  virtual void apply() override
  {
    std::cout << "Applying Dummy1" << std::endl;
  }
};

// A standard block that does specify the vector index
template<typename P, typename V, std::size_t i>
class DummyBlock2
  : public Dune::BlockLab::BlockBase<P, V, i>
{
  public:
  template<typename Context>
  DummyBlock2(Context&, const Dune::ParameterTree&)
  {
    std::cout << "Block Dummy2 with i=" << i << " created!" << std::endl;
  }

  virtual ~DummyBlock2() = default;

  virtual void apply() override
  {
    std::cout << "Applying Dummy2" << std::endl;
  }
};

// A standard block that uses two consecutive vectors
template<typename P, typename V, std::size_t i>
class DummyBlock3
  : public Dune::BlockLab::BlockBase<P, V, i>
{
  public:
  static constexpr std::size_t vectors = 2;

  template<typename Context>
  DummyBlock3(Context&, const Dune::ParameterTree&)
  {
    std::cout << "Block Dummy3 with i=" << i << " created!" << std::endl;
  }

  virtual ~DummyBlock3() = default;

  virtual void apply() override
  {
    std::cout << "Applying Dummy3" << std::endl;
  }
};

// A parent block
template<typename P, typename V>
class DummyParent
  : public Dune::BlockLab::ParentBlockBase<P, V>
{
  public:
  template<typename Context>
  DummyParent(Context& ctx, const Dune::ParameterTree& config)
    : Dune::BlockLab::ParentBlockBase<P, V>(ctx, config)
  {
    std::cout << "Block DummyParent constructed" << std::endl;
  }
};

int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // A dummy input
  Dune::ParameterTree config;
  config["solver.blocks"] = "parent, dummy1";
  config["parent.blocks"] = "dummy2, dummy3";
  config["dummy2.vector"] = "second";
  config["dummy3.vector"] = "first";

  // Construct a vector type - just to have one that has the correct Traits
  using Grid = Dune::UGGrid<2>;
  using GridProvider = Dune::BlockLab::StructuredSimplexGridProvider<Grid>;
  using VectorProvider = Dune::BlockLab::PkFemVectorProvider<GridProvider, 1>;

  auto grid = std::make_shared<GridProvider>(config);
  auto vector1 = std::make_shared<VectorProvider>(grid, "first");
  auto vector2 = std::make_shared<VectorProvider>(grid, "second");

  using UserParams = std::tuple<>;

  using Context = Dune::BlockLab::ConstructionContext<UserParams, VectorProvider, VectorProvider>;

  Context ctx(helper, config, vector1, vector2);

  ctx.registerBlock<DummyBlock1>("dummy1");
  ctx.registerBlock<DummyBlock2>("dummy2");
  ctx.registerBlock<DummyBlock3>("dummy3");
  ctx.registerBlock<DummyParent>("parent");

  auto solver = ctx.constructSolver(config.sub("solver"));
  solver->apply();

  return 0;
}
