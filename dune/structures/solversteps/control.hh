#ifndef DUNE_STRUCTURES_SOLVERSTEPS_CONTROL_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_CONTROL_HH

/** Some solversteps that implement control flow within
 *  solver steps - e.g. repeating a sequence of subblocks
 */

#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/traits.hh>


template<typename... V>
class RepeatStep
  : public StepCollectionStep<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  RepeatStep(const Dune::ParameterTree& params)
    : repeats(params.get<int>("iterations", 2))
  {}

  virtual ~RepeatStep() = default;

  virtual void apply() override
  {
    for(int i=0; i<repeats; ++i)
      for(auto step : this->steps)
	step->apply();
  }

  private:
  int repeats;
};

#endif
