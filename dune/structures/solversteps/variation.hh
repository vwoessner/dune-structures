#ifndef DUNE_STRUCTURES_SOLVERSTEPS_VARIATION_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_VARIATION_HH

#include<dune/common/parametertree.hh>
#include<dune/common/shared_ptr.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/traits.hh>

#include<memory>
#include<vector>
#include<type_traits>


template<typename... V>
class ContinuousVariationTransitionStep
  : public StepCollectionStep<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  ContinuousVariationTransitionStep(std::string name, int iterations=5, double start=0.0, double end=1.0)
    : name(name)
    , iterations(iterations)
    , start(start)
    , end(end)
  {}

  virtual ~ContinuousVariationTransitionStep() {}

  virtual void set_solver(std::shared_ptr<typename Traits::Solver> solver_) override
  {
    this->solver = solver_;
    this->solver->introduce_parameter(name, start);
    for (auto step : this->steps)
      step->set_solver(solver_);
  }

  virtual void apply() override
  {
    double val = start;
    for (int i=0; i<iterations; ++i)
    {
      val += (end - start) / iterations;
      this->solver->update_parameter(name, val);
      for (auto step : this->steps)
        step->apply();
    }
  }

  private:
  std::string name;
  int iterations;
  double start, end;
};


template<typename... V>
class DiscreteVariationTransitionStep
  : public StepCollectionStep<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  DiscreteVariationTransitionStep(std::string name, std::vector<typename Traits::Parameter> values)
    : name(name)
    , values(values)
  {}

  virtual ~DiscreteVariationTransitionStep() {}

  virtual void set_solver(std::shared_ptr<typename Traits::Solver> solver_) override
  {
    this->solver = solver_;
    this->solver->introduce_parameter(name, values[0]);
    for (auto step : this->steps)
      step->set_solver(solver_);
  }

  virtual void apply() override
  {
    for (auto val: values)
    {
      this->solver->update_parameter(name, val);
      for (auto step : this->steps)
        step->apply();
    }
  }

  private:
  std::string name;
  std::vector<typename Traits::Parameter> values;
};

#endif
