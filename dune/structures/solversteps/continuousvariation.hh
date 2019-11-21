#ifndef DUNE_STRUCTURES_SOLVERSTEPS_CONTINUOUSVARIATION_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_CONTINUOUSVARIATION_HH

#include<dune/common/shared_ptr.hh>
#include<dune/structures/solversteps/base.hh>


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

  virtual void apply(typename Base::Vector& vector, typename Base::ConstraintsContainer& cc) override
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

#endif
