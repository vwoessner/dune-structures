#ifndef DUNE_STRUCTURES_SOLVERSTEPS_INTERPOLATION_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_INTERPOLATION_HH

#include<dune/pdelab.hh>
#include<dune/pdelab/gridfunctionspace/tags.hh>
#include<dune/structures/callableadapters.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/traits.hh>
#include<dune/typetree/simpletransformationdescriptors.hh>
#include<dune/typetree/utility.hh>

#include<functional>


template<typename... V>
class InterpolationTransitionStep
  : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  using FunctionSignature = typename Traits::Range(typename Traits::Entity, typename Traits::GlobalCoordinate);

  InterpolationTransitionStep(std::function<FunctionSignature> func)
  {
    funcs.fill(func);
  }

  template<typename... FUNCS,
           typename std::enable_if<Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount == sizeof...(FUNCS), int>::type = 0>
  InterpolationTransitionStep(FUNCS... funcs) : funcs{funcs...}
  {}

  InterpolationTransitionStep(const std::array<std::function<FunctionSignature>,
                                               Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount
                                               >& funcs)
    : funcs(funcs)
  {}

  virtual ~InterpolationTransitionStep() {}

  virtual void apply(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer>) override
  {
    auto& gfs = vector->gridFunctionSpace();
    auto gf = makeGridFunctionTreeFromCallables(gfs, funcs);

    std::cout << "Interpolating into solution vector" << std::endl;
    Dune::PDELab::interpolate(gf, gfs, *vector);
  }

  private:
  // Store the lambdas
  std::array<std::function<FunctionSignature>,
             Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount
             > funcs;
};

#endif
