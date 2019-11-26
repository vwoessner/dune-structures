#ifndef DUNE_STRUCTURES_SOLVERCONSTRUCTION_HH
#define DUNE_STRUCTURES_SOLVERCONSTRUCTION_HH

#include<dune/common/parametertree.hh>
#include<dune/structures/transitionsolver.hh>
#include<dune/structures/utilities.hh>

#include<functional>
#include<map>
#include<memory>
#include<string>


template<typename>
class ConstructionContext;


template<typename Step>
std::shared_ptr<TransitionSolverStepBase<typename Step::Base::Vector>> default_constructed(const ConstructionContext<typename Step::Base::Vector>&, const Dune::ParameterTree&)
{
  return std::make_shared<Step>();
}


template<typename Vector>
class ConstructionContext
{
  public:
  using StepBase = std::shared_ptr<TransitionSolverStepBase<Vector>>;
  using RegisterFunction = std::function<StepBase(const ConstructionContext<Vector>&, const Dune::ParameterTree&)>;

  ConstructionContext(Dune::MPIHelper& helper, const Dune::ParameterTree& config)
    : helper(helper)
    , rootconfig(config)
  {
    registerStep("onetoone", default_constructed<OneToOneMappingChecker<Vector>>);
  }

  void registerStep(std::string identifier, RegisterFunction func)
  {
    mapping[identifier] = func;
  }

  StepBase construct_step(const Dune::ParameterTree& config)
  {
    auto identifier = config.get<std::string>("type");
    return mapping[identifier](*this, config);
  }

  std::unique_ptr<TransitionSolver<Vector>> construct(const Dune::ParameterTree& config)
  {
    auto solver = std::make_unique<TransitionSolver<Vector>>();

    auto stepstr = config.get<std::string>("steps");
    auto steps = str_split(stepstr);

    for (auto step: steps)
    {
      str_trim(step);
      if(config.hasSub(step))
        solver->add(construct_step(config.sub(step)));
      else if(rootconfig.hasSub(step))
        solver->add(construct_step(rootconfig.sub(step)));
      else
      {
        Dune::ParameterTree dummy;
        dummy["type"] = step;
        solver->add(construct_step(dummy));
      }
    }

    return solver;
  }

  // The reference members that might be needed for construction of steps
  Dune::MPIHelper& helper;
  const Dune::ParameterTree& rootconfig;

  private:
  // The stored mapping for each step
  std::map<std::string, RegisterFunction> mapping;
};

#endif
