#ifndef DUNE_STRUCTURES_SOLVERSTEPS_BASE_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_BASE_HH

#include<dune/structures/solversteps/traits.hh>

#include<memory>
#include<string>


template<typename... V>
class TransitionSolverStepBase
{
  public:
  using Traits = SimpleStepTraits<V...>;
  using Base = TransitionSolverStepBase<V...>;

  // The virtual interface - pretty simple
  virtual ~TransitionSolverStepBase() {}

  virtual void set_solver(std::shared_ptr<typename Traits::Solver> solver_)
  {
    solver = solver_;
  }

  virtual void pre(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc)
  {}

  virtual void apply(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc)
  {}

  virtual void post(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc)
  {}

  virtual void update_parameter(std::string name, typename Traits::Parameter param)
  {}

  protected:
  std::shared_ptr<typename Traits::Solver> solver;
};


template<typename... V>
class StepCollectionStep
  : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  StepCollectionStep()
    : steps(0)
  {}

  virtual ~StepCollectionStep() {}

  virtual void set_solver(std::shared_ptr<typename Traits::Solver> solver_) override
  {
    this->solver = solver_;
    for (auto step : steps)
      step->set_solver(solver_);
  }

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
  {
    for (auto step : steps)
      step->update_parameter(name, param);
  }

  virtual void pre(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc) override
  {
    for (auto step : steps)
      step->pre(vector, cc);
  }

  virtual void apply(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc) override
  {
    for (auto step : steps)
      step->apply(vector, cc);
  }


  virtual void post(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc) override
  {
    for (auto step : steps)
      step->post(vector, cc);
  }

  template<typename STEP>
  void add(std::shared_ptr<STEP> step)
  {
    steps.push_back(step);
    step->set_solver(this->solver);
  }

  template<typename STEP>
  void add(STEP& step)
  {
    add(Dune::stackobject_to_shared_ptr(step));
  }

  protected:
  std::vector<std::shared_ptr<TransitionSolverStepBase<V...>>> steps;
};


template<typename BaseT, typename... V>
class WrapperStep
  : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  WrapperStep(std::shared_ptr<BaseT> step)
    : step(step)
  {}

  virtual ~WrapperStep() {}

  virtual void set_solver(std::shared_ptr<typename Traits::Solver> solver_) override
  {
    this->solver = solver_;
    this->step->set_solver(solver_);
  }

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
  {
    step->update_parameter(name, param);
  }

  virtual void pre(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc) override
  {
    step->pre(vector, cc);
  }

  virtual void apply(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc) override
  {
    step->apply(vector, cc);
  }

  virtual void post(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc) override
  {
    step->post(vector, cc);
  }

  protected:
  std::shared_ptr<BaseT> step;
};

#endif
