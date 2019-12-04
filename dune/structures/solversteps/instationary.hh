#ifndef DUNE_STRUCTURES_SOLVERSTEPS_INSTATIONARY_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_INSTATIONARY_HH

#include<dune/common/parametertree.hh>
#include<dune/pdelab.hh>
#include<dune/structures/solversteps/base.hh>
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
    : theta(1.0)
  {}

  virtual ~OneStepMethodStep() {}

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    auto gfs = vector->gridFunctionSpaceStorage();
    Dune::PDELab::ISTL::BCRSMatrixBackend<> mb(21);
    sgo = std::make_shared<SpatialGridOperator>(*gfs, *cc, *gfs, *cc, *slop, mb);
    tgo = std::make_shared<TemporalGridOperator>(*gfs, *cc, *gfs, *cc, *tlop, mb);

    linearsolver = std::make_shared<LinearSolver>(0);
    newton = std::make_shared<NewtonSolver>(*igo, *vector, *linearsolver);
    newton->setVerbosityLevel(2);

    onestepmethod = std::make_shared<OneStepMethod>(theta, *igo, *newton);
  }

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
//    onestepmethod->apply(vector);
  }

  void set_spatial_localoperator(std::shared_ptr<SLOP> lop)
  {
    slop = lop;
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
  Dune::PDELab::OneStepThetaParameter<double> theta;
};


template<typename Vector>
class InstationarySolverStep
  : public VariationTransitionStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  InstationarySolverStep(const Dune::ParameterTree& config)
    : dt(config.get<double>("dt")), Tend(config.get<double>("Tend"))
  {}

  InstationarySolverStep(double dt, double Tend)
    : dt(dt), Tend(Tend)
  {}

  virtual ~InstationarySolverStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    double time = 0.0;
    while (time < Tend)
    {
      this->update_transition_value(time);

      // Apply the solver
      for (auto step : this->steps)
        step->apply(vector, cc);

      time += dt;
    }
  }

  private:
  double dt;
  double Tend;
};

#endif
