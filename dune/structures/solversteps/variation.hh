#ifndef DUNE_STRUCTURES_SOLVERSTEPS_VARIATION_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_VARIATION_HH

#include<dune/common/parametertree.hh>
#include<dune/common/shared_ptr.hh>
#include<dune/structures/solversteps/base.hh>

#include<memory>
#include<vector>
#include<type_traits>


// Forward declarations
template<typename Vector>
class ParametrizedMaterialStepBase;

template<typename Vector>
class MaterialDependantStepBase;


template<typename Vector>
class VariationTransitionStepBase
  : public StepCollectionStep<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  template<typename STEP>
  void add(std::shared_ptr<STEP> step)
  {
    if constexpr (std::is_convertible<STEP*, MaterialDependantStepBase<Vector>*>::value)
      this->steps.push_back(std::make_shared<ParametrizedMaterialStepBase<Vector>>(step, [](auto& c, auto p) { return c; } ));
    else
      this->steps.push_back(step);
  }

  template<typename STEP>
  void add(STEP& step)
  {
    add(Dune::stackobject_to_shared_ptr(step));
  }
};


template<typename Vector>
class ContinuousVariationTransitionStep
  : public VariationTransitionStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  ContinuousVariationTransitionStep(std::string name, int iterations=5, double start=0.0, double end=1.0)
    : name(name)
    , iterations(iterations)
    , start(start)
    , end(end)
  {}

  virtual ~ContinuousVariationTransitionStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    double val = start;
    for (int i=0; i<iterations; ++i)
    {
      val += (end - start) / iterations;
      this->update_parameter(name, val);
      for (auto step : this->steps)
        step->apply(vector, cc);
    }
  }

  private:
  std::string name;
  int iterations;
  double start, end;
};


template<typename Vector>
class DiscreteVariationTransitionStep
  : public VariationTransitionStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  DiscreteVariationTransitionStep(std::string name, std::vector<typename Base::Parameter> values)
    : name(name)
    , values(values)
  {}

  virtual ~DiscreteVariationTransitionStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    for (auto val: values)
    {
      this->update_parameter(name, val);
      for (auto step : this->steps)
        step->apply(vector, cc);
    }
  }

  private:
  std::string name;
  std::vector<typename Base::Parameter> values;
};

#endif
