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
    return P(tree.get<double>("value"));
  else if (type == "int")
    return P(tree.get<int>("value"));
  else if (type == "string")
    return P(tree.get<std::string>("value"));
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
                      std::shared_ptr<typename StepTraits::Grid> grid,
                      typename StepTraits::EntitySet es,
                      std::shared_ptr<std::vector<int>> physical)
    : helper(helper)
    , rootconfig(config)
    , grid(grid)
    , es(es)
    , physical(physical)
    , solver(std::make_unique<TransitionSolver<V...>>(grid, es))
  {
    // Parse
    {
      auto vecstr = config.get<std::string>("vectors", "solution");
      auto vectors  = str_split(vecstr);

      std::size_t i = 0;
      for (auto vec: vectors)
      {
        str_trim(vec);
        vector_name_to_index[vec] = i++;
      }
    }

    registerStep("adaptivity",
		 [](auto& ctx, const auto& p)
		 {
                   auto step = std::make_shared<AdaptivitySolverStep<V...>>();
                   ctx.add_children(step, p);
                   return step;
		 });

    registerVectorStep("constraints",
                 [](auto i, const auto& ctx, const auto& p)
                 {
                   return std::make_shared<ConstraintsTransitionStep<i, V...>>(
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

    registerVectorStep("elasticity_operator",
                 [](auto i, const auto& ctx, const auto& p)
                 {
                   return std::make_shared<ElasticityOperatorStep<i, V...>>(p);
                 });

    registerVectorStep("elasticity_mass_operator",
                 [](auto i, const auto& ctx, const auto& p)
                 {
                   return std::make_shared<ElasticityMassOperatorStep<i, V...>>(p);
                 });

    if constexpr (dim == 2)
      registerVectorStep("fibre_operator",
                  [](auto i, auto& ctx, const auto& p)
                  {
                    return std::make_shared<FibreReinforcedElasticitySolverStep<i, V...>>(ctx.rootconfig, p);
                  });

    registerStep("fibre_refinement",
		 [](const auto& ctx, const auto& p)
		 {
                    return std::make_shared<FiberVicinityMarkerStep<V...>>();
		 });

    registerStep("filelogger",
		 [](const auto& ctx, const auto& p)
		 {
                    return std::make_shared<FileLoggerStep<V...>>(p);
		 });

    registerVectorStep("interpolation",
                 [](auto i, const auto& ctx, const auto& p)
                 {
                   return std::make_shared<InterpolationTransitionStep<i, V...>>(
                     get_callable_array<LocalEntitySignature, V...>(*ctx.solver, p["functions"]));
                 });

    registerVectorStep("linearsolver",
                 [](auto i, const auto& ctx, const auto& p)
                 {
                   return std::make_shared<LinearSolverStep<i, V...>>(p);
                 });

    registerStep("material",
                 [](const auto& ctx, const auto& p)
                 {
                   return std::make_shared<MaterialInitialization<V...>>(ctx.es, ctx.physical, p, ctx.rootconfig);
                 });

    registerVectorStep("newton",
                 [](auto i, auto& ctx, const auto& p)
                 {
                   return std::make_shared<NewtonSolverTransitionStep<i, V...>>(p);
                 });

    registerVectorStep("onestep",
                 [](auto i, auto& ctx, const auto& p)
                 {
                   return std::make_shared<OneStepMethodStep<i, V...>>(p);
                 });

    registerVectorStep("onetoone",
                 [](auto i, const auto& ctx, const auto& p)
                 {
                   return std::make_shared<OneToOneMappingChecker<i, V...>>();
                 });

    registerStep("parameter",
                 [](const auto& ctx, const auto& p)
                 {
                   auto param = dynamic_parse<typename StepTraits::Parameter>(p);
                   return std::make_shared<ParameterSetup<V...>>(p.template get<std::string>("name"), param);
                 });

    registerVectorStep("probe",
                 [](auto i, const auto& ctx, const auto& p)
                 {
                   return std::make_shared<ProbeTransitionStep<i, V...>>(ctx.es.gridView(), p);
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
    registerVectorStep("transformation",
                 [](auto i, auto& ctx, const auto& p)
                 {
                   return std::make_shared<TransformationTransitionStep<i, V...>>(get_transformation<Coord(Coord, Coord), V...>(*ctx.solver, p.template get<std::string>("functions")));
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

    registerVectorStep("vis_solution",
                 [](auto i, const auto& ctx, const auto& p)
                 {
                   return std::make_shared<SolutionVisualizationStep<i, V...>>();
                 });

    registerVectorStep("vis_vonmises",
                 [](auto i, const auto& ctx, const auto& p)
                 {
                   return std::make_shared<VonMisesStressVisualizationStep<i, V...>>(p);
                 });

    registerStep("vis_fibredistance",
                 [](const auto& ctx, const auto& p)
                 {
                   return std::make_shared<FibreDistanceVisualizationStep<V...>>(p, ctx.rootconfig);
                 });

    registerStep("vis_indexset",
                 [](const auto& ctx, const auto& p)
				 {
                   return std::make_shared<IndexSetVisualizationStep<V...>>();
				 });
  }

  template<typename Func>
  void registerVectorStep(std::string identifier, Func&& func)
  {
    is_vector[identifier] = true;
    Dune::Hybrid::forEach(Dune::Hybrid::integralRange(std::integral_constant<std::size_t, 0>{},
                                                      std::integral_constant<std::size_t, sizeof...(V)>{}),
                          [this, identifier, func](auto i){
                            auto subident = identifier + "_" + std::to_string(i);
                            this->mapping[subident] = [i, func](auto& ctx, const auto& p){
                              return func(i, ctx, p);
                            };
                          });
  }

  template<typename Func>
  void registerStep(std::string identifier, Func&& func)
  {
    is_vector[identifier] = false;
    mapping[identifier] = std::forward<Func>(func);
  }

  StepBasePointer construct_step(std::string stepname, const Dune::ParameterTree& config)
  {
    auto identifier = config.get<std::string>("type", stepname);
    if (is_vector[identifier])
      identifier = identifier + "_" + std::to_string(vector_name_to_index[config.get<std::string>("vector", "solution")]);
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
    solver->setVectors(vectors);
    solver->setConstraintsContainers(constraints_containers);
    add_children(solver, config);
    return solver;
  }

  void setVectors(const std::shared_ptr<V>&... vectors_)
  {
    vectors = {vectors_...};
  }

  void setConstraintsContainers(const std::shared_ptr<typename VectorToConstraintsContainer<V>::type>&... containers_)
  {
    constraints_containers = {containers_...};
  }

  template<std::size_t i=0>
  auto getVector() const
  {
    return std::get<i>(vectors);
  }

  template<std::size_t i=0>
  auto getConstraintsContainer() const
  {
    return std::get<i>(constraints_containers);
  }

  // The reference members that might be needed for construction of steps
  Dune::MPIHelper& helper;
  const Dune::ParameterTree& rootconfig;
  std::shared_ptr<typename StepTraits::Grid> grid;
  typename StepTraits::EntitySet es;
  std::shared_ptr<std::vector<int>> physical;

  // The vectors and constraints containers
  std::tuple<std::shared_ptr<V>...> vectors;
  std::tuple<std::shared_ptr<typename VectorToConstraintsContainer<V>::type>...> constraints_containers;

  // The stored mapping for each step
  std::map<std::string, RegisterFunction> mapping;
  std::shared_ptr<TransitionSolver<V...>> solver;
  std::map<std::string, std::size_t> vector_name_to_index;
  std::map<std::string, bool> is_vector;
};

#endif
