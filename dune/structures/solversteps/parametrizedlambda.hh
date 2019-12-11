#ifndef DUNE_STRUCTURES_SOLVERSTEPS_PARAMETRIZEDLAMBDA_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_PARAMETRIZEDLAMBDA_HH

#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/variation.hh>

#include<functional>
#include<tuple>


template<typename Signature, typename Param>
class ParametrizedLambda
{};


template<typename Range, typename... Arguments, typename Param>
class ParametrizedLambda<Range(Arguments...), Param>
{
  public:
  template<typename FUNC>
  ParametrizedLambda(const FUNC& func)
    : param(std::make_shared<Param>())
    , func(func)
  {}

  Range operator()(const Arguments&... args) const
  {
    return func(args..., *param);
  }

  void set_parameter(Param p)
  {
    *param = p;
  }

  private:
  std::shared_ptr<Param> param;
  std::function<Range(Arguments..., Param)> func;
};


template<typename Step>
struct PFuncMixin
{
  using Base = TransitionSolverStepBase<typename Step::Base::Vector>;
  using FunctionSignature = typename Step::FunctionSignature;
  std::shared_ptr<ParametrizedLambda<FunctionSignature, typename Base::Parameter>> pfunc;
};


template<typename Step>
class ParametrizedLambdaVariationTransitionStepBase
  : public PFuncMixin<Step>,
    public WrapperStep<typename Step::Base::Vector>
{
  using PFuncMixin<Step>::pfunc;

  public:
  using Base = TransitionSolverStepBase<typename Step::Base::Vector>;
  using FunctionSignature = typename Step::FunctionSignature;

  template<typename FUNC>
  ParametrizedLambdaVariationTransitionStepBase(std::string name, const FUNC& func)
    : PFuncMixin<Step>{std::make_shared<ParametrizedLambda<FunctionSignature, typename Base::Parameter>>(func)}
    , WrapperStep<typename Base::Vector>(std::make_shared<Step>(ParametrizedLambda<FunctionSignature, typename Base::Parameter>(*pfunc)))
    , myname(name)
  {}

  virtual ~ParametrizedLambdaVariationTransitionStepBase() {}

  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    if (name == myname)
      pfunc->set_parameter(param);
  }

  private:
  std::string myname;
};

#endif
