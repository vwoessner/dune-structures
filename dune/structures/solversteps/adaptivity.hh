#ifndef DUNE_STRUCTURES_SOLVERSTEPS_ADAPTIVITY_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_ADAPTIVITY_HH

#include<dune/pdelab.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/traits.hh>


template<typename... V>
class AdaptivitySolverStep
  : public StepCollectionStep<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  virtual ~AdaptivitySolverStep() = default;

  virtual void apply() override
  {
    // TODO: Here, we should check that the grid is unmarked.
    //       All marking should be done by substeps of this one.

    // Recursively call all substeps. These are expected to mark cells for refinement.
    for (auto step : this->steps)
      step->apply();

    std::cout << "Adapting the grid..." << std::endl;
    // TODO: This currently adapts only one vector where of course it should adapt all of them
    auto vector = this->solver->getVector();
    auto gfs = vector->gridFunctionSpaceStorage();
    Dune::PDELab::adapt_grid(*(this->solver->getGrid()), *gfs, *vector, 4);
    // TODO: Integration order must come from somewhere
  }
};


template<typename... V>
class FiberVicinityMarkerStep
  : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  virtual ~FiberVicinityMarkerStep() = default;

  virtual void apply() override
  {
    std::cout << "Marking elements in the vicinity of a fibre for refinement..." << std::endl;
    auto vector = this->solver->getVector();
    auto es = vector->gridFunctionSpaceStorage()->gridView();

    // Pseudo global refine for now
    for (auto e : elements(es))
      this->solver->getGrid()->mark(1, e);

    // This is where the actual marking happens
  }
};

#endif
