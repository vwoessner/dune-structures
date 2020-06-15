#ifndef DUNE_BLOCKLAB_BLOCKS_BLOCKBASE_HH
#define DUNE_BLOCKLAB_BLOCKS_BLOCKBASE_HH

#include<dune/blocklab/blocks/blocktraits.hh>
#include<dune/blocklab/utilities/stringsplit.hh>
#include<dune/common/parametertree.hh>
#include<dune/common/shared_ptr.hh>

#include<memory>
#include<tuple>
#include<vector>


namespace Dune::BlockLab {

  // This does not have an implementation - only the specialization for std::tuple's has.
  template<typename P, typename V>
  class AbstractBlockBase;


  template<typename... P, typename... V>
  class AbstractBlockBase<std::tuple<P...>, std::tuple<V...>>
    : public std::enable_shared_from_this<AbstractBlockBase<std::tuple<P...>, std::tuple<V...>>>
  {
    public:
    using BlockBase = AbstractBlockBase<std::tuple<P...>, std::tuple<V...>>;
    using Solver = Dune::BlockLab::BlockSolver<std::tuple<P...>, std::tuple<V...>>;
    static constexpr std::size_t vectors = 1;

    // The virtual interface - pretty simple
    virtual ~AbstractBlockBase() = default;

    virtual void set_parent(std::shared_ptr<BlockBase> parent) = 0;
    virtual void set_solver(std::shared_ptr<Solver> solver) = 0;
    virtual void setup() = 0;
    virtual void apply() = 0;
    virtual void update_parameter(std::string name, typename Solver::Parameter param) = 0;

    protected:
    std::shared_ptr<Solver> solver;
    std::shared_ptr<BlockBase> parent;
  };


  template<typename P, typename V, std::size_t i=0>
  class BlockBase
    : public AbstractBlockBase<P, V>
  {
    public:
    using Traits = BlockTraits<P, V, i>;

    virtual ~BlockBase() = default;

    virtual void set_parent(std::shared_ptr<typename Traits::BlockBase> parent) override final
    {
      this->parent = parent;
    }

    virtual void set_solver(std::shared_ptr<typename Traits::Solver> solver_) override final
    {
      this->solver = solver_;
    }

    virtual void setup() override
    {}

    virtual void apply() override
    {}

    virtual void update_parameter(std::string name, typename Traits::Parameter param) override
    {}
  };


  template<typename P, typename V>
  class ParentBlockBase
    : public AbstractBlockBase<P, V>
  {
    public:
    using Traits = BlockTraits<P, V>;

    // We delete the default constructor of this class to make sure that all
    // Context/Config constructors from derived classes call the above Context/Config
    // constructor. If you manually construct a solver, pass a (potentially empty)
    // vector of blocks instead.
    ParentBlockBase() = delete;

    template<typename Context>
    ParentBlockBase(Context& ctx, const Dune::ParameterTree& config)
    {
      auto blocks = config.get<std::string>("blocks", "");
      for (auto block: string_split(blocks))
      {
        if(config.hasSub(block))
          add(ctx.constructBlock(block, config.sub(block)));
        else if(ctx.getRootConfig().hasSub(block))
          add(ctx.constructBlock(block, ctx.getRootConfig().sub(block)));
        else
        {
          Dune::ParameterTree dummy;
          dummy["type"] = block;
          add(ctx.constructBlock(block, dummy));
        }
      }
    }

    ParentBlockBase(std::vector<std::shared_ptr<AbstractBlockBase<P, V>>> blocks)
      : blocks(blocks)
    {}

    virtual ~ParentBlockBase() = default;

    virtual void set_parent(std::shared_ptr<typename Traits::BlockBase> parent) override final
    {
      this->parent = parent;
      for (auto block : blocks)
	block->set_parent(this->shared_from_this());
    }

    virtual void set_solver(std::shared_ptr<typename Traits::Solver> solver_) override final
    {
      this->solver = solver_;
      for (auto block : blocks)
        block->set_solver(solver_);
    }

    virtual void update_parameter(std::string name, typename Traits::Parameter param) override
    {
      for (auto block : blocks)
        block->update_parameter(name, param);
    }

    virtual void apply() override
    {
      for (auto block : blocks)
        block->apply();
    }


    virtual void setup() override
    {
      for (auto block : blocks)
        block->setup();
    }

    template<typename STEP>
    void add(std::shared_ptr<STEP> block)
    {
      blocks.push_back(block);
    }

    template<typename STEP>
    void add(STEP& block)
    {
      add(Dune::stackobject_to_shared_ptr(block));
    }

    protected:
    std::vector<std::shared_ptr<AbstractBlockBase<P, V>>> blocks;
  };

} // namespace Dune::BlockLab

#endif
