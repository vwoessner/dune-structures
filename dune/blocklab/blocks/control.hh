#ifndef DUNE_BLOCKLAB_BLOCKS_CONTROL_HH
#define DUNE_BLOCKLAB_BLOCKS_CONTROL_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/common/parametertree.hh>

namespace Dune::BlockLab {

  template<typename P, typename V>
  class RepeatBlock
    : public ParentBlockBase<P, V>
  {
    public:
    template<typename Context>
    RepeatBlock(Context& ctx, const Dune::ParameterTree& config)
      : ParentBlockBase<P, V>(ctx, config)
      , repeats(config.get<int>("iterations", 2))
    {}

    virtual ~RepeatBlock() = default;

    virtual void apply() override
    {
      for(int i=0; i<repeats; ++i)
        for(auto block : this->blocks)
          block->apply();
    }

    private:
    int repeats;
  };

} // namespace Dune::BlockLab

#endif
