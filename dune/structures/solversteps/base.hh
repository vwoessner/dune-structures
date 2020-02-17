#ifndef DUNE_STRUCTURES_SOLVERSTEPS_BASE_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_BASE_HH

#include<dune/structures/material.hh>

#include<memory>
#include<string>
#include<type_traits>
#include<variant>
#include<vector>


// Forward declaration of the solver class
template<typename V>
class TransitionSolver;


template<typename V>
class TransitionSolverStepBase
{
  public:
  // Export some types from the PDELab Vector type
  using Base = TransitionSolverStepBase<V>;
  using Vector = V;
  using GridFunctionSpace = typename Vector::GridFunctionSpace;
  using GridView = typename GridFunctionSpace::Traits::GridViewType;
  static constexpr int dim = GridView::dimension;
  using Grid = typename GridView::Traits::Grid;
  using EntitySet = typename GridFunctionSpace::Traits::EntitySet;
  using ctype = typename Grid::ctype;
  using Entity = typename GridView::template Codim<0>::Entity;
  using GlobalCoordinate = typename Entity::Geometry::GlobalCoordinate;
  using Range = typename Vector::field_type;
  using ConstraintsContainer = typename GridFunctionSpace::template ConstraintsContainer<Range>::Type;
  using VectorBackend = typename GridFunctionSpace::Traits::Backend;
  using Solver = TransitionSolver<V>;

  // The possible types for parametrization of solver steps
  using Material = MaterialCollection<typename Base::EntitySet, double>;
  using Parameter = std::variant<bool,
                                 double,
                                 int,
                                 std::string,
                                 std::shared_ptr<Material>,
                                 std::shared_ptr<std::vector<int>>,
                                 Dune::ParameterTree>;

  // The virtual interface - pretty simple
  virtual ~TransitionSolverStepBase() {}

  virtual void set_solver(std::shared_ptr<Solver> solver_)
  {
    solver = solver_;
  }

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<ConstraintsContainer> cc)
  {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<ConstraintsContainer> cc)
  {}

  virtual void post(std::shared_ptr<Vector> vector, std::shared_ptr<ConstraintsContainer> cc)
  {}

  virtual void update_parameter(std::string name, Parameter param)
  {}

  protected:
  std::shared_ptr<Solver> solver;
};


template<typename Vector>
class StepCollectionStep
  : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  StepCollectionStep()
    : steps(0)
  {}

  virtual ~StepCollectionStep() {}

  virtual void set_solver(std::shared_ptr<typename Base::Solver> solver_) override
  {
    this->solver = solver_;
    for (auto step : steps)
      step->set_solver(solver_);
  }

  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    for (auto step : steps)
      step->update_parameter(name, param);
  }

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    for (auto step : steps)
      step->pre(vector, cc);
  }

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    for (auto step : steps)
      step->apply(vector, cc);
  }


  virtual void post(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    for (auto step : steps)
      step->post(vector, cc);
  }

  template<typename STEP>
  void add(std::shared_ptr<STEP> step)
  {
    steps.push_back(step);
    step->set_solver(this->solver);
  }

  template<typename STEP>
  void add(STEP& step)
  {
    add(Dune::stackobject_to_shared_ptr(step));
  }

  protected:
  std::vector<std::shared_ptr<TransitionSolverStepBase<Vector>>> steps;
};


template<typename Vector, typename BaseT=TransitionSolverStepBase<Vector>>
class WrapperStep
  : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  WrapperStep(std::shared_ptr<BaseT> step)
    : step(step)
  {}

  virtual ~WrapperStep() {}

  virtual void set_solver(std::shared_ptr<typename Base::Solver> solver_) override
  {
    this->solver = solver_;
    this->step->set_solver(solver_);
  }

  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    step->update_parameter(name, param);
  }

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    step->pre(vector, cc);
  }

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    step->apply(vector, cc);
  }

  virtual void post(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    step->post(vector, cc);
  }

  protected:
  std::shared_ptr<BaseT> step;
};

#endif
