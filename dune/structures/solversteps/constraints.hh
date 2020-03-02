#ifndef DUNE_STRUCTURES_SOLVERSTEPS_CONSTRAINTS_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_CONSTRAINTS_HH


#include<dune/pdelab.hh>
#include<dune/structures/callableadapters.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/traits.hh>

#include<functional>


template<std::size_t i, typename... V>
class ConstraintsTransitionStep
  : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = VectorStepTraits<i, V...>;
  using FunctionSignature = bool(typename Traits::GridView::Intersection, typename Traits::GridView::Intersection::Geometry::LocalCoordinate);

  ConstraintsTransitionStep(std::function<FunctionSignature> func)
  {
    funcs.fill(func);
  }

  template<typename... FUNCS,
           typename std::enable_if<Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount == sizeof...(FUNCS), int>::type = 0>
  ConstraintsTransitionStep(FUNCS... funcs) : funcs{funcs...}
  {}

  ConstraintsTransitionStep(const std::array<std::function<FunctionSignature>,
                                             Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount
                                             >& funcs)
    : funcs(funcs)
  {}

  virtual ~ConstraintsTransitionStep() {}

  virtual void apply() override
  {
    auto vector = this->solver->template getVector<i>();
    auto constraintscontainer = this->solver->template getConstraintsContainer<i>();
    auto& gfs = vector->gridFunctionSpace();
    auto bctype = makeBoundaryConditionTreeFromCallables(gfs, funcs);
    Dune::PDELab::constraints(bctype, gfs, *constraintscontainer);
    std::cout << "Assembled constraints - " << constraintscontainer->size() << " of " << gfs.size() << " dofs constrained!" << std::endl;
  }

  private:
  // Store the lambdas
  std::array<std::function<FunctionSignature>,
             Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount
             > funcs;
};

#endif
