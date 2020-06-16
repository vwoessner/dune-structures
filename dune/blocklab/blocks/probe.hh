#ifndef DUNE_BLOCKLAB_BLOCKS_PROBE_HH
#define DUNE_BLOCKLAB_BLOCKS_PROBE_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/common/parametertree.hh>
#include<dune/pdelab/common/functionutilities.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>

#include<memory>
#include<string>


namespace Dune::BlockLab {

  template<typename P, typename V, std::size_t i>
  class ProbeBlock
    : public BlockBase<P, V, i>
  {
    public:
    using Traits = BlockTraits<P, V, i>;

    using DGF = Dune::PDELab::DiscreteGridFunction<typename Traits::GridFunctionSpace, typename Traits::Vector>;
    using Probe = Dune::PDELab::GridFunctionProbe<DGF>;

    template<typename Context>
    ProbeBlock(Context&, const Dune::ParameterTree& config)
      : ProbeBlock(config.get<typename Traits::GlobalCoordinate>("position"),
		   config.get<std::string>("name", ""))
    {}

    ProbeBlock(typename Traits::GlobalCoordinate position,
	       std::string name = "")
      : name(name)
      , position(position)
    {}

    virtual ~ProbeBlock() = default;

    virtual void setup() override
    {
      if (name != "")
        this->solver->introduce_parameter(name, typename DGF::Traits::RangeType(0.0));

      auto vector = this->solver->template getVector<i>();
      probe = std::make_shared<Probe>(vector->gridFunctionSpace().gridView(), position);
    }

    virtual void apply() override
    {
      auto vector = this->solver->template getVector<i>();
      DGF dgf(vector->gridFunctionSpaceStorage(), vector);
      probe->setGridFunction(dgf);

      typename DGF::Traits::RangeType eval(0.0);
      probe->eval(eval);

      // Report the value of the probe into the solvers parameter system
      if (name != "")
        this->solver->update_parameter(name, eval);

      std::cout << "Probe " << name << ": " << eval << std::endl;
    }

    private:
    std::string name;
    typename Traits::GlobalCoordinate position;
    std::shared_ptr<Probe> probe;
  };

} // namespace Dune::BlockLab

#endif
