#ifndef DUNE_STRUCTURES_SOLVERSTEPS_DISCRETEVARIATION_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_DISCRETEVARIATION_HH

#include<dune/structures/solversteps/base.hh>


template<typename Vector, typename ValueType>
class DiscreteVariationTransitionStep
  : public VariationTransitionStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  DiscreteVariationTransitionStep(std::vector<ValueType> values)
    : VariationTransitionStepBase<Vector>(), values(values)
  {}

  virtual ~DiscreteVariationTransitionStep() {}

  virtual void apply(typename Base::Vector& vector, typename Base::ConstraintsContainer& cc)
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
