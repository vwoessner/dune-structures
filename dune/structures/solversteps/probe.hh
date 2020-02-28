#ifndef DUNE_STRUCTURES_SOLVERSTEPS_PROBE_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_PROBE_HH

#include<dune/common/parametertree.hh>
#include<dune/pdelab.hh>
#include<dune/structures/solversteps/base.hh>

#include<memory>


template<typename... V>
class ProbeTransitionStep
  : public TransitionSolverStepBase<V...>
{
  public:
  using Base = TransitionSolverStepBase<V...>;
  using DGF = Dune::PDELab::VectorDiscreteGridFunction<typename Base::GridFunctionSpace, typename Base::Vector>;
  using Probe = Dune::PDELab::GridFunctionProbe<DGF>;

  ProbeTransitionStep(const typename Base::GridView& gv, const Dune::ParameterTree& config)
    : config(config)
    , probe(gv, config.get<typename Base::GlobalCoordinate>("position"))
  {}

  virtual ~ProbeTransitionStep() {}

  virtual void apply(std::shared_ptr<typename Base::Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    DGF dgf(vector->gridFunctionSpaceStorage(), vector);
    probe.setGridFunction(dgf);

    typename DGF::Traits::DomainType eval(0.0);
    probe.eval(eval);

    if (config.get<bool>("displace", false))
      eval += config.get<typename Base::GlobalCoordinate>("position");

    std::cout << "Probe " << config.get<std::string>("name", "") << ": " << eval << std::endl;
  }

  private:
  Dune::ParameterTree config;
  Probe probe;
};

#endif
