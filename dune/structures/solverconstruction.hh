#ifndef DUNE_STRUCTURES_SOLVERCONSTRUCTION_HH
#define DUNE_STRUCTURES_SOLVERCONSTRUCTION_HH

#include<dune/common/parametertree.hh>
#include<dune/structures/muparser.hh>
#include<dune/structures/transitionsolver.hh>
#include<dune/structures/utilities.hh>
#include<dune/typetree/utility.hh>

#include<array>
#include<functional>
#include<map>
#include<memory>
#include<string>


template<typename>
class ConstructionContext;


template<typename Step>
std::shared_ptr<TransitionSolverStepBase<typename Step::Base::Vector>> default_constructed(const ConstructionContext<typename Step::Base::Vector>&, const Dune::ParameterTree&)
{
  return std::make_shared<Step>();
}


template<typename Signature>
std::function<Signature> get_callable(std::string expr)
{
  return MuParserCallable<Signature>(expr);
}


template<typename GFS, typename Signature>
std::array<std::function<Signature>,
           Dune::TypeTree::TreeInfo<GFS>::leafCount> get_callable_array(std::string expr)
{
  constexpr auto len = Dune::TypeTree::TreeInfo<GFS>::leafCount;
  std::array<std::function<Signature>, len> result;

  auto exprs = str_split(expr);
  if (exprs.size() == 1)
    result.fill(get_callable<Signature>(exprs[0]));
  else
   std::transform(exprs.begin(), exprs.end(), result.begin(), get_callable<Signature>);

  return std::move(result);
}



template<typename Vector>
class ConstructionContext
{
  public:
  using StepBase = TransitionSolverStepBase<Vector>;
  using StepBasePointer = std::shared_ptr<StepBase>;
  using RegisterFunction = std::function<StepBase(const ConstructionContext<Vector>&, const Dune::ParameterTree&)>;

  ConstructionContext(Dune::MPIHelper& helper,
                      const Dune::ParameterTree& config,
                      typename StepBase::EntitySet es,
                      std::shared_ptr<std::vector<int>> physical)
    : helper(helper)
    , rootconfig(config)
    , es(es)
    , physical(physical)
  {
    registerStep("constraints",
                 [](const auto& ctx, auto& p)
                 {
                   return std::make_shared<ConstraintsTransitionStep<Vector>>(
                     get_callable_array<typename Vector::GridFunctionSpace, bool(Dune::FieldVector<double, 3>)>(p["functions"]));
                 });

    registerStep("elasticity",
                 [](const auto& ctx, auto& p)
                 {
                   return std::make_shared<ElasticitySolverStep<Vector>>(
                     ctx.es, ctx.physical, ctx.rootparams);
                 });

    registerStep("interpolation",
                 [](const auto& ctx, auto& p)
                 {
                   return std::make_shared<InterpolationTransitionStep<Vector>>(
                     get_callable_array<typename Vector::GridFunctionSpace, double(Dune::FieldVector<double, 3>)>(p["functions"]));
                 });

    registerStep("onetoone",
                 default_constructed<OneToOneMappingChecker<Vector>>);

    registerStep("visualization",
                 default_constructed<VisualizationStep<Vector, false>>);

  }

  void registerStep(std::string identifier, RegisterFunction func)
  {
    mapping[identifier] = func;
  }

  StepBasePointer construct_step(std::string stepname, const Dune::ParameterTree& config)
  {
    auto identifier = config.get<std::string>("type", stepname);
    return mapping[identifier](*this, config);
  }

  std::unique_ptr<TransitionSolver<Vector>> construct(const Dune::ParameterTree& config)
  {
    auto solver = std::make_unique<TransitionSolver<Vector>>();

    auto stepstr = config.get<std::string>("steps", "");

    if (stepstr == "")
      return solver;

    auto steps = str_split(stepstr);

    for (auto step: steps)
    {
      std::cout << "Trying to build a step " << step << std::endl;
      str_trim(step);
      if(config.hasSub(step))
        solver->add(construct_step(step, config.sub(step)));
      else if(rootconfig.hasSub(step))
        solver->add(construct_step(step, rootconfig.sub(step)));
      else
      {
        Dune::ParameterTree dummy;
        dummy["type"] = step;
        solver->add(construct_step(step, dummy));
      }
    }

    return solver;
  }

  // The reference members that might be needed for construction of steps
  Dune::MPIHelper& helper;
  const Dune::ParameterTree& rootconfig;
  typename StepBase::EntitySet es;
  std::shared_ptr<std::vector<int>> physical;

  private:
  // The stored mapping for each step
  std::map<std::string, RegisterFunction> mapping;
};

#endif
