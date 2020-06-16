#ifndef DUNE_BLOCKLAB_BLOCKS_VARIATION_HH
#define DUNE_BLOCKLAB_BLOCKS_VARIATION_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/construction/dataparser.hh>
#include<dune/common/parametertree.hh>

#include<string>


namespace Dune::BlockLab {

  template<typename P, typename V>
  class ContinuousVariationBlock
    : public ParentBlockBase<P, V>
  {
    public:
    using Traits = BlockTraits<P, V>;

    template<typename Context>
    ContinuousVariationBlock(Context& ctx, const Dune::ParameterTree& config)
      : ParentBlockBase<P, V>(ctx, config)
      , name(config.get<std::string>("name"))
      , iterations(config.get<int>("iterations"))
      , start(config.get<double>("start", 0.0))
      , end(config.get<double>("end", 1.0))
    {}

    virtual ~ContinuousVariationBlock() = default;

    virtual void setup() override
    {
      this->solver->introduce_parameter(name, start);

      for (auto block : this->blocks)
        block->setup();
    }

    virtual void apply() override
    {
      double val = start;
      for (int i=0; i<iterations; ++i)
      {
        val += (end - start) / iterations;
        this->solver->update_parameter(name, val);
        for (auto block : this->blocks)
          block->apply();
      }
    }

    private:
    std::string name;
    int iterations;
    double start, end;
  };


  template<typename P, typename V>
  class DiscreteVariationBlock
    : public ParentBlockBase<P, V>
  {
    public:
    using Traits = BlockTraits<P, V>;

    template<typename Context>
    DiscreteVariationBlock(Context& ctx, const Dune::ParameterTree& config)
      : ParentBlockBase<P,V>(ctx, config)
      , name(config.get<std::string>("name"))
      , values(parse_parameter_list<typename Traits::Parameter>(config))
    {}

    virtual ~DiscreteVariationBlock() = default;

    virtual void setup() override
    {
      this->solver->introduce_parameter(name, values[0]);

      for (auto block : this->blocks)
        block->setup();
    }

    virtual void apply() override
    {
      for (auto val: values)
      {
        this->solver->update_parameter(name, val);
        for (auto block : this->blocks)
          block->apply();
      }
    }

    private:
    std::string name;
    std::vector<typename Traits::Parameter> values;
  };

} // namespace Dune::BlockLab

#endif
