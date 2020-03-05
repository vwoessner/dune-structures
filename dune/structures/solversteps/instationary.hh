#ifndef DUNE_STRUCTURES_SOLVERSTEPS_INSTATIONARY_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_INSTATIONARY_HH

#include<dune/common/parametertree.hh>
#include<dune/pdelab.hh>
#include<dune/structures/callableadapters.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/interpolation.hh>
#include<dune/structures/solversteps/newton.hh>
#include<dune/structures/solversteps/traits.hh>
#include<dune/structures/solversteps/variation.hh>

#include<memory>


template<std::size_t i, typename... V>
class OneStepMethodStep
  : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = VectorStepTraits<i, V...>;

  using VirtLocalOperator = AbstractLocalOperatorInterface<typename Traits::GridFunctionSpace>;


  using GridOperator = Dune::PDELab::GridOperator<typename Traits::GridFunctionSpace,
                                                  typename Traits::GridFunctionSpace,
                                                  VirtLocalOperator,
                                                  Dune::PDELab::ISTL::BCRSMatrixBackend<>,
                                                  typename Traits::ctype,
                                                  typename Traits::Range,
                                                  typename Traits::Range,
                                                  typename Traits::ConstraintsContainer,
                                                  typename Traits::ConstraintsContainer>;

  using InstationaryGridOperator = Dune::PDELab::OneStepGridOperator<GridOperator, GridOperator>;

  using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_UMFPack;
  using NewtonSolver = Dune::PDELab::Newton<InstationaryGridOperator, LinearSolver, typename Traits::Vector>;
  using OneStepMethod = Dune::PDELab::OneStepMethod<double, InstationaryGridOperator, NewtonSolver, typename Traits::Vector>;

  OneStepMethodStep(const Dune::ParameterTree& params)
    : params(params)
    , theta(std::make_shared<Dune::PDELab::OneStepThetaParameter<double>>(params.get<double>("theta", 1.0)))
    , time(0.0)
    , timestep(0.0)
  {}

  virtual ~OneStepMethodStep() {}

  virtual void apply() override
  {
    std::cout << "Applying the one step method" << std::endl;
    auto vector = this->solver->template getVector<i>();
    onestepmethod->apply(time, timestep, *vector, *swapvector);
    vector.swap(swapvector);
  }

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
  {
    if (name == "time")
      time = std::get<double>(param);
    if (name == "timestep")
      timestep = std::get<double>(param);
    if ((name == params.get<std::string>("spatial_operator")) || (name == params.get<std::string>("temporal_operator")))
    {
      // One of the operator changed - rebuild the one step method object
      auto vector = this->solver->template getVector<i>();
      auto cc = this->solver->template getConstraintsContainer<i>();
      auto gfs = vector->gridFunctionSpaceStorage();
      Dune::PDELab::ISTL::BCRSMatrixBackend<> mb(21);
      auto slop = this->solver->template param<std::shared_ptr<VirtLocalOperator>>(params.template get<std::string>("spatial_operator"));
      auto tlop = this->solver->template param<std::shared_ptr<VirtLocalOperator>>(params.template get<std::string>("temporal_operator"));
      auto sgo = std::make_shared<GridOperator>(*gfs, *cc, *gfs, *cc, *slop, mb);
      auto tgo = std::make_shared<GridOperator>(*gfs, *cc, *gfs, *cc, *tlop, mb);
      igo = std::make_shared<InstationaryGridOperator>(*sgo, *tgo);

      linearsolver = std::make_shared<LinearSolver>(0);
      newton = std::make_shared<NewtonSolver>(*igo, *vector, *linearsolver);
      newton->setVerbosityLevel(2);
      swapvector = std::make_shared<typename Traits::Vector>(*vector);
      onestepmethod = std::make_shared<OneStepMethod>(*theta, *igo, *newton);
    }
  }

  protected:
  std::shared_ptr<GridOperator> sgo;
  std::shared_ptr<GridOperator> tgo;
  std::shared_ptr<InstationaryGridOperator> igo;
  std::shared_ptr<LinearSolver> linearsolver;
  std::shared_ptr<NewtonSolver> newton;
  std::shared_ptr<OneStepMethod> onestepmethod;

  Dune::ParameterTree params;
  std::shared_ptr<Dune::PDELab::OneStepThetaParameter<double>> theta;
  std::shared_ptr<typename Traits::Vector> swapvector;
  double time;
  double timestep;
};


template<std::size_t i, typename... V>
class VariableBoundaryOneStepMethodStep
  : public OneStepMethodStep<i, V...>
{
  public:
  using Traits = VectorStepTraits<i, V...>;

  using FunctionSignature = typename Traits::Range(typename Traits::GlobalCoordinate);

  VariableBoundaryOneStepMethodStep(const Dune::ParameterTree& params, std::function<FunctionSignature> singlefunc)
    : OneStepMethodStep<i, V...>(params)
  {
    funcs.fill(singlefunc);
  }

  template<typename... FUNCS,
           typename std::enable_if<Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount == sizeof...(FUNCS), int>::type = 0>
  VariableBoundaryOneStepMethodStep(const Dune::ParameterTree& params, FUNCS... funcs)
    : OneStepMethodStep<i, V...>(params)
    , funcs{funcs...}
  {}

  VariableBoundaryOneStepMethodStep(const Dune::ParameterTree& params,
                                    const std::array<std::function<FunctionSignature>,
                                                     Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount>& funcs)
    : OneStepMethodStep<i, V...>(params)
    , funcs(funcs)
  {}

  virtual ~VariableBoundaryOneStepMethodStep() {}

  virtual void apply() override
  {
    auto vector = this->solver->template getVector<i>();
    auto& gfs = vector->gridFunctionSpace();
    auto func = makeInstationaryGridFunctionTreeFromCallables(*this->solver, gfs, funcs);

    this->onestepmethod->apply(this->time, this->timestep, *vector, func, *this->swapvector);
    vector.swap(this->swapvector);
  }

  private:
  // Store the lambdas
  std::array<std::function<FunctionSignature>,
             Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount
             > funcs;
};


template<typename... V>
class InstationarySolverStep
  : public ContinuousVariationTransitionStep<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  InstationarySolverStep(double Tstart, double Tend, double dt)
    : ContinuousVariationTransitionStep<V...>("time")
    , dt(dt), Tstart(Tstart), Tend(Tend)
  {}

  virtual ~InstationarySolverStep() {}

  virtual void set_solver(std::shared_ptr<typename Traits::Solver> solver_)
  {
    this->solver = solver_;
    this->solver->introduce_parameter("time", Tstart);
    this->solver->introduce_parameter("timestep", dt);
  }

  virtual void apply() override
  {
    std::cout << "Starting time stepping loop" << std::endl;
    double time = Tstart;
    while (time < Tend - 1e-8)
    {
      this->update_parameter("time", time);
      this->update_parameter("timestep", dt);
      std::cout << "Performing time step " << time << " -> " << time + dt  << " with " << this->steps.size() << " steps" << std::endl;

      // Apply the solver
      for (auto step : this->steps)
        step->apply();

      time += dt;
    }
  }

  private:
  double dt;
  double Tstart;
  double Tend;
};

#endif
