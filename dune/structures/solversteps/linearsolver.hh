#ifndef DUNE_STRUCTURES_SOLVERSTEPS_LINEARSOLVER_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_LINEARSOLVER_HH

#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/traits.hh>


template<std::size_t i, typename... V>
class LinearSolverStep
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

  using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_UMFPack;
  using StationaryLinearProblemSolver = Dune::PDELab::StationaryLinearProblemSolver<GridOperator, LinearSolver, typename Traits::Vector>;

  LinearSolverStep(const Dune::ParameterTree& params)
    : params(params)
  {}

  virtual ~LinearSolverStep() {}

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
  {
    if (name == params.get<std::string>("operator"))
    {
      // The operator changed - rebuild the Newton solver object!
      auto vector = this->solver->template getVector<i>();
      auto cc = this->solver->template getConstraintsContainer<i>();
      auto gfs = vector->gridFunctionSpaceStorage();
      auto localoperator = this->solver->template param<std::shared_ptr<VirtLocalOperator>>(params.get<std::string>("operator"));
      Dune::PDELab::ISTL::BCRSMatrixBackend<> mb(21);
      gridoperator = std::make_shared<GridOperator>(*gfs, *cc, *gfs, *cc, *localoperator, mb);
      linearsolver = std::make_shared<LinearSolver>(0);
      slp = std::make_shared<StationaryLinearProblemSolver>(*gridoperator, *linearsolver, *vector, 1e-12);
    }
  }

  virtual void apply() override
  {
    std::cout << "Applying linear Solver!" << std::endl;
    slp->apply();
  }

  protected:
  Dune::ParameterTree params;
  std::shared_ptr<LinearSolver> linearsolver;
  std::shared_ptr<GridOperator> gridoperator;
  std::shared_ptr<StationaryLinearProblemSolver> slp;
};
#endif
