#ifndef DUNE_BLOCKLAB_BLOCKS_INSTATIONARY_HH
#define DUNE_BLOCKLAB_BLOCKS_INSTATIONARY_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/common/parametertree.hh>

#include<memory>

namespace Dune::BlockLab {

  template<typename P, typename V>
  class TimestepperBlock
    : public ParentBlockBase<P, V>
  {
    public:
    using Traits = BlockTraits<P, V>;

    template<typename Context>
    TimestepperBlock(Context& ctx, const Dune::ParameterTree& config)
      : ParentBlockBase<P, V>(ctx, config)
      , dt(config.get<double>("timestep"))
      , Tstart(config.get<double>("starttime"))
      , Tend(config.get<double>("endtime"))
    {}

    virtual ~TimestepperBlock() = default;

    virtual void setup() override
    {
      this->solver->introduce_parameter("time", Tstart);
      this->solver->introduce_parameter("timestep", dt);
    }

    virtual void apply() override
    {
      std::cout << "Starting time stepping loop" << std::endl;
      double time = Tstart;
      while (time < Tend - 1e-8)
      {
        this->update_parameter("time", time);
        this->update_parameter("timestep", dt);
        std::cout << "Performing time step " << time << " -> " << time + dt  << " with " << this->steps.size() << " steps" << std::endl;

        // Apply the solver
        for (auto step : this->steps)
          step->apply();

        time += dt;
      }
    }

    private:
    double dt;
    double Tstart;
    double Tend;
  };

} // namespace Dune::BlockLab

#endif
