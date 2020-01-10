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
#include<vector>


template<typename>
class ConstructionContext;


template<typename Step>
std::shared_ptr<TransitionSolverStepBase<typename Step::Base::Vector>> default_constructed(const ConstructionContext<typename Step::Base::Vector>&, const Dune::ParameterTree&)
{
  return std::make_shared<Step>();
}


template<typename Step>
std::shared_ptr<TransitionSolverStepBase<typename Step::Base::Vector>> with_tree(const ConstructionContext<typename Step::Base::Vector>&, const Dune::ParameterTree& tree)
{
  return std::make_shared<Step>(tree);
}


template<typename Vector, typename Signature>
std::function<Signature> get_callable(TransitionSolver<Vector>& solver, std::string expr)
{
  return MuParserCallable<Signature>(solver, expr);
}

template<typename Vector, typename Signature>
std::function<Signature> get_transformation(TransitionSolver<Vector>& solver, std::string expr)
{
  return MuParserTransformation<Signature>(solver, expr);
}


template<typename Vector, typename Signature>
std::array<std::function<Signature>,
           Dune::TypeTree::TreeInfo<typename Vector::GridFunctionSpace>::leafCount> get_callable_array(TransitionSolver<Vector>& solver, std::string expr)
{
  using GFS = typename Vector::GridFunctionSpace;
  constexpr auto len = Dune::TypeTree::TreeInfo<GFS>::leafCount;
  std::array<std::function<Signature>, len> result;

  auto exprs = str_split(expr);
  if (exprs.size() == 1)
    result.fill(get_callable<Vector, Signature>(solver, exprs[0]));
  else
    std::transform(exprs.begin(), exprs.end(), result.begin(), [&solver](auto it){ return get_callable<Vector, Signature>(solver, it); });

  return std::move(result);
}


template<typename P>
P dynamic_parse(const Dune::ParameterTree& tree, std::string default_type = "double")
{
  auto type = tree.get("datatype", default_type);

  if (type == "double")
    return tree.get<double>("value");
  else if (type == "int")
    return tree.get<int>("value");
  else if (type == "string")
    return tree.get<std::string>("value");
  else
    DUNE_THROW(Dune::Exception, "Cannot parse parameter");
}

template<typename P>
std::vector<P> dynamic_list_parse(const Dune::ParameterTree& tree, std::string default_type = "string")
{
  auto type = tree.get("datatype", default_type);
  auto valuestr = tree.get<std::string>("values");
  auto values = str_split(valuestr);
  std::vector<P> ret;

  for (auto val : values)
  {
    Dune::ParameterTree p;
    p["datatype"] = type;
    p["value"] = val;
    ret.push_back(dynamic_parse<P>(p));
  }

  return ret;
}

template<typename Vector>
class ConstructionContext
{
  public:
  using StepBase = TransitionSolverStepBase<Vector>;
  using StepBasePointer = std::shared_ptr<StepBase>;
  using RegisterFunction = std::function<StepBasePointer(ConstructionContext<Vector>&, const Dune::ParameterTree&)>;
  using Coord = typename StepBase::GlobalCoordinate;

  ConstructionContext(Dune::MPIHelper& helper,
                      const Dune::ParameterTree& config,
                      typename StepBase::EntitySet es,
                      std::shared_ptr<std::vector<int>> physical)
    : helper(helper)
    , rootconfig(config)
    , es(es)
    , physical(physical)
    , solver(std::make_unique<TransitionSolver<Vector>>())
  {
    registerStep("constraints",
                 [](const auto& ctx, const auto& p)
                 {
                   return std::make_shared<ConstraintsTransitionStep<Vector>>(
                     get_callable_array<Vector, bool(Coord)>(*ctx.solver, p.template get<std::string>("functions")));
                 });

    registerStep("continuousvariation",
                 [](auto& ctx, const auto& p)
                 {
                   auto step = std::make_shared<ContinuousVariationTransitionStep<Vector>>(p.template get<std::string>("parameter"),
                                                                                           p.template get<int>("iterations"),
                                                                                           p.template get<double>("start"),
                                                                                           p.template get<double>("end"));
                   step->set_solver(ctx.solver);
                   ctx.add_children(step, p);
                   return step;
                 });

    registerStep("discretevariation",
                 [](auto& ctx, const auto& p)
                 {
                   auto values = dynamic_list_parse<typename StepBase::Parameter>(p);
                   auto step = std::make_shared<DiscreteVariationTransitionStep<Vector>>(p.template get<std::string>("name"), values);
                   step->set_solver(ctx.solver);
                   ctx.add_children(step, p);
                   return step;
                 });

    registerStep("discretematerialvariation",
                 [](auto& ctx, const auto& p)
                 {
                   auto values = dynamic_list_parse<typename StepBase::Parameter>(p);
                   auto key = p.template get<std::string>("key");
                   auto name = p.template get<std::string>("name");
                   auto step = std::make_shared<DiscreteMaterialVariationTransitionStep<Vector>>(name,
                                                                                                 [key](auto& tree, auto param){ tree[key] = std::get<std::string>(param); },
                                                                                                 values);
                   step->set_solver(ctx.solver);
                   ctx.add_children(step, p);
                   return step;
                 });

    registerStep("interpolation",
                 [](const auto& ctx, const auto& p)
                 {
                   return std::make_shared<InterpolationTransitionStep<Vector>>(
                     get_callable_array<Vector, double(Coord)>(*ctx.solver, p["functions"]));
                 });

    registerStep("material",
                 [](const auto& ctx, const auto& p)
                 {
                   return std::make_shared<MaterialInitialization<Vector>>(ctx.es, ctx.physical, p, ctx.rootconfig);
                 });

    registerStep("parameter",
                 [](const auto& ctx, const auto& p)
                 {
                   auto param = dynamic_parse<typename StepBase::Parameter>(p);
                   return std::make_shared<ParameterSetup<Vector>>(p.template get<std::string>("name"), param);
                 });

    registerStep("timeloop",
                 [](auto& ctx, const auto& p)
                 {
                   auto step = std::make_shared<InstationarySolverStep<Vector>>(p.template get<double>("Tstart"),
                                                                                p.template get<double>("Tend"),
                                                                                p.template get<double>("dt"));
                   step->set_solver(ctx.solver);
                   ctx.add_children(step, p);
                   return step;
                 });
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

  template<typename ROOT>
  void add_children(std::shared_ptr<ROOT> root, const Dune::ParameterTree& config)
  {
    auto stepstr = config.get<std::string>("steps", "");

    if (stepstr == "")
      return;

    auto steps = str_split(stepstr);

    for (auto step: steps)
    {
      str_trim(step);
      if(config.hasSub(step))
        root->add(construct_step(step, config.sub(step)));
      else if(rootconfig.hasSub(step))
        root->add(construct_step(step, rootconfig.sub(step)));
      else
      {
        Dune::ParameterTree dummy;
        dummy["type"] = step;
        root->add(construct_step(step, dummy));
      }
    }
  }

  std::shared_ptr<TransitionSolver<Vector>> construct(const Dune::ParameterTree& config)
  {
    add_children(solver, config);
    return solver;
  }

  // The reference members that might be needed for construction of steps
  Dune::MPIHelper& helper;
  const Dune::ParameterTree& rootconfig;
  typename StepBase::EntitySet es;
  std::shared_ptr<std::vector<int>> physical;

  // The stored mapping for each step
  std::map<std::string, RegisterFunction> mapping;
  std::shared_ptr<TransitionSolver<Vector>> solver;
};



template<typename Vector>
class ElasticityConstructionContext
  : public ConstructionContext<Vector>
{
  public:
  using Base = ConstructionContext<Vector>;

  template<typename... Params>
  ElasticityConstructionContext(Params&&... params)
    : ConstructionContext<Vector>(std::forward<Params>(params)...)
  {
    this->registerStep("elasticity",
                       with_tree<ElasticitySolverStep<Vector>>);

    this->registerStep("onetoone",
                       default_constructed<OneToOneMappingChecker<Vector>>);

    this->registerStep("probe",
                       [](const auto& ctx, const auto& p)
                       {
                         return std::make_shared<ProbeTransitionStep<Vector>>(ctx.es.gridView(), p);
                       });

    this->registerStep("transformation",
                 [](auto& ctx, const auto& p)
                 {
                   return std::make_shared<TransformationTransitionStep<Vector>>(get_transformation<Vector, typename Base::Coord(typename Base::Coord, typename Base::Coord)>(*ctx.solver, p.template get<std::string>("functions")));
                 });

    this->registerStep("visualization",
                       [](auto& ctx, const auto& p)
                       {
                         auto step = std::make_shared<VisualizationStep<Vector, false>>(p);
                         ctx.add_children(step, p);
                         return step;
                       });

    this->registerStep("vis_mpirank",
                       [](auto& ctx, const auto& p)
                       {
                         return std::make_shared<MPIRankVisualizationStep<Vector, true>>(ctx.helper);
                       });

    this->registerStep("vis_physicalentity",
                       [](auto& ctx, const auto& p)
                       {
                         return std::make_shared<PhysicalEntityVisualizationStep<Vector, false>>(ctx.physical);
                       });

    this->registerStep("vis_solution",
                       default_constructed<SolutionVisualizationStep<Vector, false>>);

    this->registerStep("vis_vonmises",
                       default_constructed<VonMisesStressVisualizationStep<Vector, false>>);

    this->registerStep("vis_fibredistance",
                       [](const auto& ctx, const auto& p)
                       {
                         return std::make_shared<FibreDistanceVisualizationStep<Vector>>(p, ctx.rootconfig);
                       });
  }
};

template<typename Vector>
class ElastodynamicsConstructionContext
    : public ConstructionContext<Vector>
  {
    public:
    using Base = ConstructionContext<Vector>;

    template<typename... Params>
    ElastodynamicsConstructionContext(Params&&... params)
      : ConstructionContext<Vector>(std::forward<Params>(params)...)
    {
      this->registerStep("elastodynamics",
                         [](const auto& ctx, const auto& p)
                         {
                           if (p.hasKey("functions"))
                             return std::make_shared<ElastoDynamicsSolverStep<Vector>>(p, get_callable_array<Vector, double(typename Base::Coord)>(*ctx.solver, p["functions"]));
                           else
                             return std::make_shared<ElastoDynamicsSolverStep<Vector>>(p);
                         });

      this->registerStep("visualization",
                         [](auto& ctx, const auto& p)
                         {
                           auto step = std::make_shared<VisualizationStep<Vector, true>>(p);
                           ctx.add_children(step, p);
                           return step;
                         });

      this->registerStep("vis_mpirank",
                         [](auto& ctx, const auto& p)
                         {
                           return std::make_shared<MPIRankVisualizationStep<Vector, true>>(ctx.helper);
                         });

      this->registerStep("vis_physicalentity",
                         [](auto& ctx, const auto& p)
                         {
                           return std::make_shared<PhysicalEntityVisualizationStep<Vector, true>>(ctx.physical);
                         });

      this->registerStep("vis_solution",
                         default_constructed<SolutionVisualizationStep<Vector, true>>);
    }
  };

#endif
