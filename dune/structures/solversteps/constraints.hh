#ifndef DUNE_STRUCTURES_SOLVERSTEPS_CONSTRAINTS_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_CONSTRAINTS_HH


#include<dune/pdelab.hh>
#include<dune/structures/callableadapters.hh>
#include<dune/structures/solversteps/base.hh>

#include<functional>


template<typename Vector>
class ConstraintsTransitionStep : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;
  using FunctionSignature = bool(typename Base::GlobalCoordinate);

  ConstraintsTransitionStep(std::function<FunctionSignature> func)
  {
    funcs.fill(func);
  }

  template<typename... FUNCS,
           typename std::enable_if<Dune::TypeTree::TreeInfo<typename Base::GridFunctionSpace>::leafCount == sizeof...(FUNCS), int>::type = 0>
  ConstraintsTransitionStep(FUNCS... funcs) : funcs{funcs...}
  {}

  ConstraintsTransitionStep(const std::array<std::function<FunctionSignature>,
                                             Dune::TypeTree::TreeInfo<typename Base::GridFunctionSpace>::leafCount
                                             >& funcs)
    : funcs(funcs)
  {}

  virtual ~ConstraintsTransitionStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> constraintscontainer) override
  {
    auto& gfs = vector->gridFunctionSpace();
    auto bctype = makeBoundaryConditionTreeFromCallables(gfs, funcs);
    Dune::PDELab::constraints(bctype, gfs, *constraintscontainer);
    std::cout << "Assembled constraints - " << constraintscontainer->size() << " of " << gfs.size() << " dofs constrained!" << std::endl;
  }

  private:
  // Store the lambdas
  std::array<std::function<FunctionSignature>,
             Dune::TypeTree::TreeInfo<typename Base::GridFunctionSpace>::leafCount
             > funcs;
};

#endif
