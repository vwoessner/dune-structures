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

    /** TODO:
     * In the below implementation, only one vector is treated, where all vectors from
     * the variadic template list V... should be considered. Unfortunately all these must
     * be handled in one call of adapt_grid. The code should look something like this
     * in Python generator syntax, but I yet need to understand how this looks like in C++
     * metaprogramming world.
     *
     * adapt_grid(grid,
     *            transferSolutions(*this->solver->template getVector<i>()->gridFunctionSpaceStorage(),
     *                              *this->solver->template getVector<i>(),
     *                              4)
     *            for i in range(sizeof...(V)));
     *
     * Of course all the below hassle with const_cast was removed from this minimal example.
     *
     * Another thing that bothers me, where I do not understand PDELab semantics enough:
     * If two of the vectors in V... use the same grid function space, is it okay to pass
     * that gfs twice into adapt or do I actually need to group my vectors by gfs object???
     * I doubt I will ever succeed in doing that in meta-programming, I would much rather
     * duplicate and invade the PDELab adaptivity interface.
     */

    auto vector = this->solver->getVector();
    auto gfs = vector->gridFunctionSpaceStorage();

    // The extraction of the grid function space is super-ugly here: PDELab requires a non-const
    // reference to the grid function space, but our way of not redundantly carrying around gfs's
    // gives us only shared pointers to `const GFS`. We solve this by casting away the const.
    Dune::PDELab::adapt_grid(*(this->solver->getGrid()),
			     *const_cast<typename std::remove_const<typename decltype(gfs)::element_type>::type*>(gfs.get()),
			     *vector,
			     4);

    // TODO: Integration order was hardcoded to 4 above, but must come from somewhere
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
