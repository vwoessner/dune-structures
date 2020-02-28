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


template<typename... V>
class ConstructionContext;


template<typename Step, typename... V>
std::shared_ptr<TransitionSolverStepBase<V...>> default_constructed(const ConstructionContext<V...>&, const Dune::ParameterTree&)
{
  return std::make_shared<Step>();
}


template<typename Step, typename... V>
std::shared_ptr<TransitionSolverStepBase<V...>> with_tree(const ConstructionContext<V...>&, const Dune::ParameterTree& tree)
{
  return std::make_shared<Step>(tree);
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


template<typename... V>
class ConstructionContext
{
  public:
  using StepBase = TransitionSolverStepBase<V...>;
  using StepBasePointer = std::shared_ptr<StepBase>;
  using RegisterFunction = std::function<StepBasePointer(ConstructionContext<V...>&, const Dune::ParameterTree&)>;
  using StepTraits = SimpleStepTraits<V...>;
  using Entity = typename StepTraits::Entity;
  using Coord = typename StepTraits::GlobalCoordinate;
  using LocalCoord = typename Entity::Geometry::LocalCoordinate;
  using Intersection = typename StepTraits::GridView::Intersection;
  using IntersectionLocalCoord = typename Intersection::Geometry::LocalCoordinate;
  static constexpr int dim = StepTraits::dim;

  using LocalEntitySignature = double(Entity, LocalCoord);
  using LocalIntersectionSignature = bool(Intersection, IntersectionLocalCoord);

  ConstructionContext(Dune::MPIHelper& helper,
                      const Dune::ParameterTree& config,
                      typename StepTraits::EntitySet es,
                      std::shared_ptr<std::vector<int>> physical)
    : helper(helper)
    , rootconfig(config)
    , es(es)
    , physical(physical)
    , solver(std::make_unique<TransitionSolver<V...>>(es))
  {
    registerStep("constraints",
                 [](const auto& ctx, const auto& p)
                 {
                   return std::make_shared<ConstraintsTransitionStep<V...>>(
                     get_callable_array<LocalIntersectionSignature, V...>(*ctx.solver, p.template get<std::string>("functions")));
                 });

    registerStep("continuousvariation",
                 [](auto& ctx, const auto& p)
                 {
                   auto step = std::make_shared<ContinuousVariationTransitionStep<V...>>(p.template get<std::string>("name"),
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
                   auto values = dynamic_list_parse<typename StepTraits::Parameter>(p);
                   auto step = std::make_shared<DiscreteVariationTransitionStep<V...>>(p.template get<std::string>("name"), values);
                   step->set_solver(ctx.solver);
                   ctx.add_children(step, p);
                   return step;
                 });

    registerStep("discretematerialvariation",
                 [](auto& ctx, const auto& p)
                 {
                   auto values = dynamic_list_parse<typename StepTraits::Parameter>(p);
                   auto key = p.template get<std::string>("key");
                   auto name = p.template get<std::string>("name");
                   auto step = std::make_shared<DiscreteMaterialVariationTransitionStep<V...>>(name,
                                                                                                 [key](auto& tree, auto param){ tree[key] = std::get<std::string>(param); },
                                                                                                 values);
                   step->set_solver(ctx.solver);
                   ctx.add_children(step, p);
                   return step;
                 });

    registerStep("elasticity",
                 with_tree<ElasticitySolverStep<V...>, V...>);

    if constexpr (dim == 2)
      registerStep("fibrereinforcedelasticity",
                   [](auto& ctx, const auto& p)
                   {
                     return std::make_shared<FibreReinforcedElasticitySolverStep<0, V...>>(ctx.rootconfig, p);
                   });

    registerStep("interpolation",
                 [](const auto& ctx, const auto& p)
                 {
                   return std::make_shared<InterpolationTransitionStep<V...>>(
                     get_callable_array<LocalEntitySignature, V...>(*ctx.solver, p["functions"]));
                 });

    registerStep("material",
                 [](const auto& ctx, const auto& p)
                 {
                   return std::make_shared<MaterialInitialization<V...>>(ctx.es, ctx.physical, p, ctx.rootconfig);
                 });

    registerStep("onetoone",
                 default_constructed<OneToOneMappingChecker<V...>, V...>);

    registerStep("parameter",
                 [](const auto& ctx, const auto& p)
                 {
                   auto param = dynamic_parse<typename StepTraits::Parameter>(p);
                   return std::make_shared<ParameterSetup<V...>>(p.template get<std::string>("name"), param);
                 });

    registerStep("probe",
                 [](const auto& ctx, const auto& p)
                 {
                   return std::make_shared<ProbeTransitionStep<V...>>(ctx.es.gridView(), p);
                 });

    registerStep("quasistatic_elasticity",
                 [](const auto& ctx, const auto& p)
                 {
                   return std::make_shared<QuasiStaticElastoDynamicsSolverStep<V...>>(ctx.rootconfig);
                 });

    registerStep("timeloop",
                 [](auto& ctx, const auto& p)
                 {
                   auto step = std::make_shared<InstationarySolverStep<V...>>(p.template get<double>("Tstart"),
                                                                                p.template get<double>("Tend"),
                                                                                p.template get<double>("dt"));
                   step->set_solver(ctx.solver);
                   ctx.add_children(step, p);
                   return step;
                 });

    if constexpr (false) {
    registerStep("transformation",
                 [](auto& ctx, const auto& p)
                 {
                   return std::make_shared<TransformationTransitionStep<V...>>(get_transformation<Coord(Coord, Coord), V...>(*ctx.solver, p.template get<std::string>("functions")));
                 });
    }

    registerStep("visualization",
                 [](auto& ctx, const auto& p)
                 {
                   auto step = std::make_shared<VisualizationStep<V...>>(p);
                   ctx.add_children(step, p);
                   return step;
                 });

    registerStep("vis_mpirank",
                 [](auto& ctx, const auto& p)
                 {
                   return std::make_shared<MPIRankVisualizationStep<V...>>(ctx.helper);
                 });

    registerStep("vis_physicalentity",
                 [](auto& ctx, const auto& p)
                 {
                   return std::make_shared<PhysicalEntityVisualizationStep<V...>>(ctx.physical);
                 });

    registerStep("vis_solution",
                 default_constructed<SolutionVisualizationStep<V...>, V...>);

    registerStep("vis_vonmises",
                 default_constructed<VonMisesStressVisualizationStep<V...>, V...>);

    registerStep("vis_fibredistance",
                 [](const auto& ctx, const auto& p)
                 {
                   return std::make_shared<FibreDistanceVisualizationStep<V...>>(p, ctx.rootconfig);
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

  std::shared_ptr<TransitionSolver<V...>> construct(const Dune::ParameterTree& config)
  {
    add_children(solver, config);
    return solver;
  }

  // The reference members that might be needed for construction of steps
  Dune::MPIHelper& helper;
  const Dune::ParameterTree& rootconfig;
  typename StepTraits::EntitySet es;
  std::shared_ptr<std::vector<int>> physical;

  // The stored mapping for each step
  std::map<std::string, RegisterFunction> mapping;
  std::shared_ptr<TransitionSolver<V...>> solver;
};

#endif
