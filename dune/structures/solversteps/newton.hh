#ifndef DUNE_STRUCTURES_SOLVERSTEPS_NEWTON_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_NEWTON_HH

#include<dune/pdelab.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/traits.hh>


template<typename LocalOperator, typename... V>
class NewtonSolverTransitionStep
  : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  using GridOperator = Dune::PDELab::GridOperator<typename Traits::GridFunctionSpace,
                                                  typename Traits::GridFunctionSpace,
                                                  LocalOperator,
                                                  Dune::PDELab::ISTL::BCRSMatrixBackend<>,
                                                  typename Traits::ctype,
                                                  typename Traits::Range,
                                                  typename Traits::Range,
                                                  typename Traits::ConstraintsContainer,
                                                  typename Traits::ConstraintsContainer>;

  using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_UMFPack;
  using NewtonSolver = Dune::PDELab::Newton<GridOperator, LinearSolver, typename Traits::Vector>;

  NewtonSolverTransitionStep(const Dune::ParameterTree& params)
    : params(params)
    , localoperator(0)
  {}

  NewtonSolverTransitionStep(const Dune::ParameterTree& params, LocalOperator& localoperator)
    : params(params)
    , localoperator(Dune::stackobject_to_shared_ptr(localoperator))
  {}

  NewtonSolverTransitionStep(const Dune::ParameterTree& params, std::shared_ptr<LocalOperator> localoperator)
    : params(params)
    , localoperator(localoperator)
  {}

  virtual ~NewtonSolverTransitionStep() {}

  void set_localoperator(std::shared_ptr<LocalOperator> lop)
  {
    localoperator = lop;
  }

  std::shared_ptr<LocalOperator> get_localoperator()
  {
    return localoperator;
  }

  virtual void pre(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc) override
  {
    std::cout << "Building Newton" << std::endl;
    auto gfs = vector->gridFunctionSpaceStorage();
    Dune::PDELab::ISTL::BCRSMatrixBackend<> mb(21);
    gridoperator = std::make_shared<GridOperator>(*gfs, *cc, *gfs, *cc, *localoperator, mb);
    linearsolver = std::make_shared<LinearSolver>(0);
    newton = std::make_shared<NewtonSolver>(*gridoperator, *vector, *linearsolver);
    newton->setParameters(params);
    newton->setVerbosityLevel(2);
  }

  virtual void apply(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc) override
  {
    std::cout << "Applying Newton Solver!" << std::endl;
    newton->apply();
  }

  protected:
  Dune::ParameterTree params;
  std::shared_ptr<LocalOperator> localoperator;
  std::shared_ptr<LinearSolver> linearsolver;
  std::shared_ptr<GridOperator> gridoperator;
  std::shared_ptr<NewtonSolver> newton;
};

#endif
