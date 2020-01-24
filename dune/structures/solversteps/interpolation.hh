#ifndef DUNE_STRUCTURES_SOLVERSTEPS_INTERPOLATION_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_INTERPOLATION_HH

#include<dune/pdelab.hh>
#include<dune/pdelab/gridfunctionspace/tags.hh>
#include<dune/structures/callableadapters.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/typetree/simpletransformationdescriptors.hh>
#include<dune/typetree/utility.hh>

#include<functional>


template<typename Vector>
class InterpolationTransitionStep
  : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;
  using FunctionSignature = typename Base::Range(typename Base::Entity, typename Base::GlobalCoordinate);

  InterpolationTransitionStep(std::function<FunctionSignature> func)
  {
    funcs.fill(func);
  }

  template<typename... FUNCS,
           typename std::enable_if<Dune::TypeTree::TreeInfo<typename Base::GridFunctionSpace>::leafCount == sizeof...(FUNCS), int>::type = 0>
  InterpolationTransitionStep(FUNCS... funcs) : funcs{funcs...}
  {}

  InterpolationTransitionStep(const std::array<std::function<FunctionSignature>,
                                               Dune::TypeTree::TreeInfo<typename Base::GridFunctionSpace>::leafCount
                                               >& funcs)
    : funcs(funcs)
  {}

  virtual ~InterpolationTransitionStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    auto& gfs = vector->gridFunctionSpace();
    auto gf = makeGridFunctionTreeFromCallables(gfs, funcs);

    std::cout << "Interpolating into solution vector" << std::endl;
    Dune::PDELab::interpolate(gf, gfs, *vector);
  }

  private:
  // Store the lambdas
  std::array<std::function<FunctionSignature>,
             Dune::TypeTree::TreeInfo<typename Base::GridFunctionSpace>::leafCount
             > funcs;
};

#endif
