#ifndef DUNE_BLOCKLAB_BLOCKS_PARAMETER_HH
#define DUNE_BLOCKLAB_BLOCKS_PARAMETER_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/construction/dataparser.hh>
#include<dune/common/parametertree.hh>

#include<memory>

namespace Dune::BlockLab {

  template<typename P, typename V>
  class ParameterBlock
    : public BlockBase<P, V>
  {
    public:
    using Traits = BlockTraits<P, V>;

    template<typename Context>
    ParameterBlock(Context&, const Dune::ParameterTree& config)
      : ParameterBlock(config)
    {}

    ParameterBlock(const Dune::ParameterTree& config)
      : name(config.get<std::string>("name"))
      , param(parse_parameter<typename Traits::Parameter>(config))
    {}

    virtual ~ParameterBlock() = default;

    virtual void setup() override
    {
      this->solver->introduce_parameter(name, param);
    }

    private:
    std::string name;
    typename Traits::Parameter param;
  };

} // namespace Dune::BlockLab

#endif
