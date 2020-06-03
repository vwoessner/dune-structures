#ifndef DUNE_STRUCTURES_SOLVERSTEPS_ADAPTIVITY_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_ADAPTIVITY_HH

#include<dune/pdelab.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/traits.hh>

#include<tuple>


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

    /**
     * The following code is horribly complicated. So it might deserve some remarks:
     *
     * * This is necessary because all vectors to be adapted need to be passed to the
     *   `adaptGrid` function in one go.
     * *  PDELab requires a non-const reference to the grid function space, but our way
     *    of not redundantly carrying around gfs's gives us only shared pointers to `const GFS`.
     *    We solve this by casting away the const...
     * *  Another thing that bothers me, where I do not understand PDELab semantics enough:
     *    If two of the vectors in V... use the same grid function space, is it okay to pass
     *    that gfs twice into adapt or do I actually need to group my vectors by gfs object???
     *    I doubt I will ever succeed in doing that in meta-programming, I would much rather
     *    duplicate and invade the PDELab adaptivity interface.
     * *  The `4` below hardcodes the integration order. It should of course be configurable in
     *    some way, but I will not target that before I actually think about coarsening.
     */
    auto transfer = std::apply(
      [](auto... v)
      {
        return std::make_tuple(
          Dune::PDELab::transferSolutions(
            *const_cast<typename std::remove_const<typename decltype(v)::element_type::GridFunctionSpace>::type*>(v->gridFunctionSpaceStorage().get()),
            4,
            *v
          )...
        );
       },
       this->solver->getVectors()
     );

    std::apply(
      [this](auto... v)
      {
        Dune::PDELab::adaptGrid(
          *(this->solver->getGrid()),
	  v...
	);
      },
      transfer
    );
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
