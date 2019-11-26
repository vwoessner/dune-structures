#ifndef DUNE_STRUCTURES_SOLVERSTEPS_PARAMETRIZEDLAMBDA_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_PARAMETRIZEDLAMBDA_HH

#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/variation.hh>

#include<functional>
#include<tuple>


template<typename Signature, typename... Params>
class ParametrizedLambda
{};


template<typename Range, typename... Arguments, typename... Params>
class ParametrizedLambda<Range(Arguments...), Params...>
{
  public:

  template<typename FUNC>
  ParametrizedLambda(const FUNC& func)
    : params(std::make_shared<std::tuple<Params...>>()), func(func)
  {}

  Range operator()(const Arguments&... args) const
  {
    return std::apply(func, std::tuple_cat(std::make_tuple(args...), *params));
  }

  void set_parameters(const Params&... p)
  {
    *params = std::make_tuple(p...);
  }

  template<std::size_t pos>
  void set_parameter(typename std::tuple_element<pos, std::tuple<Params...>>::type p)
  {
    std::get<pos>(*params) = p;
  }

  private:
  std::shared_ptr<std::tuple<Params...>> params;
  std::function<Range(Arguments..., Params...)> func;
};


template<typename WrappedStep, typename... Params>
class ParametrizedLambdaVariationTransitionStepBase
  : public ParametrizedTransitionStepBase<typename WrappedStep::Base::Vector, Params...>
{
  public:
  using Base = TransitionSolverStepBase<typename WrappedStep::Base::Vector>;
  using FunctionSignature = typename WrappedStep::FunctionSignature;

  template<typename FUNC>
  ParametrizedLambdaVariationTransitionStepBase(FUNC func)
    : pfunc(func),
      step(std::make_shared<WrappedStep>(pfunc))
  {}

  virtual ~ParametrizedLambdaVariationTransitionStepBase() {}

  virtual void update_transition_value(Params... params) override
  {
    pfunc.set_parameters(params...);
  }

  virtual void apply(std::shared_ptr<typename Base::Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    step->apply(vector, cc);
  }

  private:
  ParametrizedLambda<FunctionSignature, Params...> pfunc;
  std::shared_ptr<TransitionSolverStepBase<typename Base::Vector>> step;
};

#endif
