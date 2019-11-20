#ifndef DUNE_STRUCTURES_SOLVERSTEPS_NEWTON_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_NEWTON_HH

#include<dune/pdelab.hh>
#include<dune/structures/onetoone.hh>
#include<dune/structures/solversteps/base.hh>


template<typename Vector, typename LocalOperator>
class NewtonSolverTransitionStep : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  using GridOperator = Dune::PDELab::GridOperator<typename Base::GridFunctionSpace,
                                                  typename Base::GridFunctionSpace,
                                                  LocalOperator,
                                                  Dune::PDELab::ISTL::BCRSMatrixBackend<>,
                                                  typename Base::ctype,
                                                  typename Base::Range,
                                                  typename Base::Range,
                                                  typename Base::ConstraintsContainer,
                                                  typename Base::ConstraintsContainer>;

  using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_UMFPack;

  NewtonSolverTransitionStep(LocalOperator& localoperator)
    : localoperator(localoperator)
  {}

  void apply(Vector& vector, typename Base::ConstraintsContainer& constraintscontainer) override
  {
    // Extract stuff
    auto gfs = vector.gridFunctionSpace();

    // Build a grid operator
    Dune::PDELab::ISTL::BCRSMatrixBackend<> mb(21);
    GridOperator gridoperator(gfs, constraintscontainer, gfs, constraintscontainer, localoperator, mb);

    // Build the Newton solver
    LinearSolver linearsolver;
    using Newton = Dune::PDELab::Newton<GridOperator, LinearSolver, Vector>;
    Newton newton(gridoperator, vector, linearsolver);

    std::cout << "Applying Newton Solver!" << std::endl;
    newton.apply();
  }

  private:
  LocalOperator& localoperator;
};

#endif
