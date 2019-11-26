#ifndef DUNE_STRUCTURES_SOLVERSTEPS_VARIATION_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_VARIATION_HH

#include<dune/common/parametertree.hh>
#include<dune/common/shared_ptr.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/construction.hh>

#include<memory>
#include<vector>
#include<type_traits>


template<typename Vector, typename... Params>
class NoopParametrizationWrapper
  : public ParametrizedTransitionStepBase<Vector, Params...>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  NoopParametrizationWrapper(std::shared_ptr<TransitionSolverStepBase<Vector>> step)
    : step(step)
  {}

  NoopParametrizationWrapper(const Dune::ParameterTree& params)
    : step(construct_step<Vector>(params))
  {}

  virtual ~NoopParametrizationWrapper() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    step->apply(vector, cc);
  }

  private:
  std::shared_ptr<TransitionSolverStepBase<Vector>> step;
};


template<typename Vector, typename... Params>
class VariationTransitionStepBase
  : public ParametrizedTransitionStepBase<Vector, Params...>
{
  public:
  using Base = ParametrizedTransitionStepBase<Vector, Params...>;

  VariationTransitionStepBase()
    : steps(0)
  {}

  virtual ~VariationTransitionStepBase() {}

  virtual void update_transition_value(Params... val) override
  {
    for (auto step : steps)
      step->update_transition_value(val...);
  }

  template<typename STEP>
  void add(std::shared_ptr<STEP> step)
  {
    if constexpr (std::is_convertible<STEP*, ParametrizedTransitionStepBase<Vector, Params...>*>::value)
      steps.push_back(step);
    else
      steps.push_back(std::make_shared<NoopParametrizationWrapper<Vector, Params...>>(step));
  }

  template<typename STEP>
  void add(STEP& step)
  {
    add(Dune::stackobject_to_shared_ptr(step));
  }

  protected:
  std::vector<std::shared_ptr<ParametrizedTransitionStepBase<Vector, Params...>>> steps;
};


template<typename Vector>
class ContinuousVariationTransitionStep
  : public VariationTransitionStepBase<Vector, double>
{
  public:
  using Base = VariationTransitionStepBase<Vector, double>;

  ContinuousVariationTransitionStep(int iterations=5, double start=0.0, double end=1.0)
    : VariationTransitionStepBase<Vector, double>(),
      iterations(iterations),
      start(start),
      end(end)
  {}

  virtual ~ContinuousVariationTransitionStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    double val = start;
    for (int i=0; i<iterations; ++i)
    {
      val += (end - start) / iterations;
      this->update_transition_value(val);
      for (auto step : this->steps)
        step->apply(vector, cc);
    }
  }

  private:
  int iterations;
  double start, end;
};


template<typename Vector, typename ValueType>
class DiscreteVariationTransitionStep
  : public VariationTransitionStepBase<Vector, ValueType>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  DiscreteVariationTransitionStep(std::vector<ValueType> values)
    : VariationTransitionStepBase<Vector, ValueType>(), values(values)
  {}

  virtual ~DiscreteVariationTransitionStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    for (auto val: values)
    {
      this->update_transition_value(val);
      for (auto step : this->steps)
        step->apply(vector, cc);
    }
  }

  private:
  std::vector<ValueType> values;
};

#endif
