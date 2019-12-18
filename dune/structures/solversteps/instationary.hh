#ifndef DUNE_STRUCTURES_SOLVERSTEPS_INSTATIONARY_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_INSTATIONARY_HH

#include<dune/common/parametertree.hh>
#include<dune/pdelab.hh>
#include<dune/structures/callableadapters.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/interpolation.hh>
#include<dune/structures/solversteps/newton.hh>
#include<dune/structures/solversteps/variation.hh>

#include<memory>


template<typename Vector, typename SLOP, typename TLOP>
class OneStepMethodStep
  : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  using SpatialGridOperator = Dune::PDELab::GridOperator<typename Base::GridFunctionSpace,
                                                         typename Base::GridFunctionSpace,
                                                         SLOP,
                                                         Dune::PDELab::ISTL::BCRSMatrixBackend<>,
                                                         typename Base::ctype,
                                                         typename Base::Range,
                                                         typename Base::Range,
                                                         typename Base::ConstraintsContainer,
                                                         typename Base::ConstraintsContainer>;

  using TemporalGridOperator = Dune::PDELab::GridOperator<typename Base::GridFunctionSpace,
                                                          typename Base::GridFunctionSpace,
                                                          TLOP,
                                                          Dune::PDELab::ISTL::BCRSMatrixBackend<>,
                                                          typename Base::ctype,
                                                          typename Base::Range,
                                                          typename Base::Range,
                                                          typename Base::ConstraintsContainer,
                                                          typename Base::ConstraintsContainer>;

  using InstationaryGridOperator = Dune::PDELab::OneStepGridOperator<SpatialGridOperator, TemporalGridOperator>;
  using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_UMFPack;
  using NewtonSolver = Dune::PDELab::Newton<InstationaryGridOperator, LinearSolver, Vector>;
  using OneStepMethod = Dune::PDELab::OneStepMethod<double, InstationaryGridOperator, NewtonSolver, Vector>;

  OneStepMethodStep()
    : theta(std::make_shared<Dune::PDELab::OneStepThetaParameter<double>>(1.0))
    , time(0.0)
    , timestep(0.0)
  {}

  virtual ~OneStepMethodStep() {}

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    auto gfs = vector->gridFunctionSpaceStorage();
    Dune::PDELab::ISTL::BCRSMatrixBackend<> mb(21);
    sgo = std::make_shared<SpatialGridOperator>(*gfs, *cc, *gfs, *cc, *slop, mb);
    tgo = std::make_shared<TemporalGridOperator>(*gfs, *cc, *gfs, *cc, *tlop, mb);
    igo = std::make_shared<InstationaryGridOperator>(*sgo, *tgo);

    linearsolver = std::make_shared<LinearSolver>(0);
    newton = std::make_shared<NewtonSolver>(*igo, *vector, *linearsolver);
    newton->setVerbosityLevel(2);

    swapvector = std::make_shared<Vector>(*vector);

    onestepmethod = std::make_shared<OneStepMethod>(*theta, *igo, *newton);
  }

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    std::cout << "Applying the one step method" << std::endl;
    onestepmethod->apply(time, timestep, *vector, *swapvector);
    vector.swap(swapvector);
  }

  virtual void update_parameter(std::string name, typename Base::Parameter param) override
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
  std::shared_ptr<Vector> swapvector;
  double time;
  double timestep;
};


template<typename Vector, typename SLOP, typename TLOP>
class VariableBoundaryOneStepMethodStep
  : public OneStepMethodStep<Vector, SLOP, TLOP>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;
  using FunctionSignature = typename Base::Range(typename Base::GlobalCoordinate);

  VariableBoundaryOneStepMethodStep(std::function<FunctionSignature> singlefunc)
  {
    funcs.fill(singlefunc);
  }

  template<typename... FUNCS,
           typename std::enable_if<Dune::TypeTree::TreeInfo<typename Base::GridFunctionSpace>::leafCount == sizeof...(FUNCS), int>::type = 0>
  VariableBoundaryOneStepMethodStep(FUNCS... funcs)
    : funcs{funcs...}
  {}

  VariableBoundaryOneStepMethodStep(const std::array<std::function<FunctionSignature>,
                                                     Dune::TypeTree::TreeInfo<typename Base::GridFunctionSpace>::leafCount>& funcs)
    : funcs(funcs)
  {}

  virtual ~VariableBoundaryOneStepMethodStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    auto& gfs = vector->gridFunctionSpace();
    auto func = makeInstationaryGridFunctionTreeFromCallables(*this->solver, gfs, funcs);

    this->onestepmethod->apply(this->time, this->timestep, *vector, func, *this->swapvector);
    vector.swap(this->swapvector);
  }

  private:
  // Store the lambdas
  std::array<std::function<FunctionSignature>,
             Dune::TypeTree::TreeInfo<typename Base::GridFunctionSpace>::leafCount
             > funcs;
};


template<typename Vector>
class InstationarySolverStep
  : public ContinuousVariationTransitionStep<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  InstationarySolverStep(double Tstart, double Tend, double dt)
    : ContinuousVariationTransitionStep<Vector>("time")
    , dt(dt), Tstart(Tstart), Tend(Tend)
  {}

  virtual ~InstationarySolverStep() {}

  virtual void set_solver(std::shared_ptr<typename Base::Solver> solver_)
  {
    this->solver = solver_;
    this->solver->introduce_parameter("time", Tstart);
    this->solver->introduce_parameter("timestep", dt);
  }

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
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
        step->apply(vector, cc);

      time += dt;
    }
  }

  private:
  double dt;
  double Tstart;
  double Tend;
};

#endif
