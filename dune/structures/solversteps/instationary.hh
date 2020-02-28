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


template<typename SLOP, typename TLOP, typename... V>
class OneStepMethodStep
  : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  using SpatialGridOperator = Dune::PDELab::GridOperator<typename Traits::GridFunctionSpace,
                                                         typename Traits::GridFunctionSpace,
                                                         SLOP,
                                                         Dune::PDELab::ISTL::BCRSMatrixBackend<>,
                                                         typename Traits::ctype,
                                                         typename Traits::Range,
                                                         typename Traits::Range,
                                                         typename Traits::ConstraintsContainer,
                                                         typename Traits::ConstraintsContainer>;

  using TemporalGridOperator = Dune::PDELab::GridOperator<typename Traits::GridFunctionSpace,
                                                          typename Traits::GridFunctionSpace,
                                                          TLOP,
                                                          Dune::PDELab::ISTL::BCRSMatrixBackend<>,
                                                          typename Traits::ctype,
                                                          typename Traits::Range,
                                                          typename Traits::Range,
                                                          typename Traits::ConstraintsContainer,
                                                          typename Traits::ConstraintsContainer>;

  using InstationaryGridOperator = Dune::PDELab::OneStepGridOperator<SpatialGridOperator, TemporalGridOperator>;
  using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_UMFPack;
  using NewtonSolver = Dune::PDELab::Newton<InstationaryGridOperator, LinearSolver, typename Traits::Vector>;
  using OneStepMethod = Dune::PDELab::OneStepMethod<double, InstationaryGridOperator, NewtonSolver, typename Traits::Vector>;

  OneStepMethodStep()
    : theta(std::make_shared<Dune::PDELab::OneStepThetaParameter<double>>(1.0))
    , time(0.0)
    , timestep(0.0)
  {}

  virtual ~OneStepMethodStep() {}

  virtual void pre() override
  {
    auto vector = this->solver->getVector();
    auto cc = this->solver->getConstraintsContainer();
    auto gfs = vector->gridFunctionSpaceStorage();
    Dune::PDELab::ISTL::BCRSMatrixBackend<> mb(21);
    sgo = std::make_shared<SpatialGridOperator>(*gfs, *cc, *gfs, *cc, *slop, mb);
    tgo = std::make_shared<TemporalGridOperator>(*gfs, *cc, *gfs, *cc, *tlop, mb);
    igo = std::make_shared<InstationaryGridOperator>(*sgo, *tgo);

    linearsolver = std::make_shared<LinearSolver>(0);
    newton = std::make_shared<NewtonSolver>(*igo, *vector, *linearsolver);
    newton->setVerbosityLevel(2);

    swapvector = std::make_shared<typename Traits::Vector>(*vector);

    onestepmethod = std::make_shared<OneStepMethod>(*theta, *igo, *newton);
  }

  virtual void apply() override
  {
    std::cout << "Applying the one step method" << std::endl;
    auto vector = this->solver->getVector();
    onestepmethod->apply(time, timestep, *vector, *swapvector);
    vector.swap(swapvector);
  }

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
  {
    if (name == "time")
      time = std::get<double>(param);
    if (name == "timestep")
      timestep = std::get<double>(param);
  }

  std::shared_ptr<SLOP> get_spatial_localoperator()
  {
    return slop;
  }

  void set_spatial_localoperator(std::shared_ptr<SLOP> lop)
  {
    slop = lop;
  }

  std::shared_ptr<TLOP> get_temporal_localoperator()
  {
    return tlop;
  }

  void set_temporal_localoperator(std::shared_ptr<TLOP> lop)
  {
    tlop = lop;
  }

  protected:
  std::shared_ptr<SLOP> slop;
  std::shared_ptr<TLOP> tlop;
  std::shared_ptr<SpatialGridOperator> sgo;
  std::shared_ptr<TemporalGridOperator> tgo;
  std::shared_ptr<InstationaryGridOperator> igo;
  std::shared_ptr<LinearSolver> linearsolver;
  std::shared_ptr<NewtonSolver> newton;
  std::shared_ptr<OneStepMethod> onestepmethod;
  std::shared_ptr<Dune::PDELab::OneStepThetaParameter<double>> theta;
  std::shared_ptr<typename Traits::Vector> swapvector;
  double time;
  double timestep;
};


template<typename SLOP, typename TLOP, typename... V>
class VariableBoundaryOneStepMethodStep
  : public OneStepMethodStep<SLOP, TLOP, V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  using FunctionSignature = typename Traits::Range(typename Traits::GlobalCoordinate);

  VariableBoundaryOneStepMethodStep(std::function<FunctionSignature> singlefunc)
  {
    funcs.fill(singlefunc);
  }

  template<typename... FUNCS,
           typename std::enable_if<Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount == sizeof...(FUNCS), int>::type = 0>
  VariableBoundaryOneStepMethodStep(FUNCS... funcs)
    : funcs{funcs...}
  {}

  VariableBoundaryOneStepMethodStep(const std::array<std::function<FunctionSignature>,
                                                     Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount>& funcs)
    : funcs(funcs)
  {}

  virtual ~VariableBoundaryOneStepMethodStep() {}

  virtual void apply() override
  {
    auto vector = this->solver->getVector();
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
