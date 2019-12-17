#ifndef DUNE_STRUCTURES_SOLVERSTEPS_PARAMETER_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_PARAMETER_HH

#include<dune/structures/solversteps/base.hh>


template<typename Vector>
class ParameterSetup
  : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  ParameterSetup(std::string name, typename Base::Parameter param)
    : name(name)
    , param(param)
  {}

  virtual ~ParameterSetup() {}

  virtual void set_solver(std::shared_ptr<typename Base::Solver> solver_) override
  {
    this->solver = solver_;
    this->solver->introduce_parameter(name, param);
  }

  private:
  std::string name;
  typename Base::Parameter param;
};

#endif
