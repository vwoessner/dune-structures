#ifndef DUNE_STRUCTURES_SOLVERSTEPS_PARAMETER_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_PARAMETER_HH

#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/traits.hh>


template<typename... V>
class ParameterSetup
  : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  ParameterSetup(std::string name, typename Traits::Parameter param)
    : name(name)
    , param(param)
  {}

  virtual ~ParameterSetup() {}

  virtual void set_solver(std::shared_ptr<typename Traits::Solver> solver_) override
  {
    this->solver = solver_;
    this->solver->introduce_parameter(name, param);
  }

  private:
  std::string name;
  typename Traits::Parameter param;
};

#endif
