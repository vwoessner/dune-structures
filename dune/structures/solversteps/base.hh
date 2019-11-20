#ifndef DUNE_STRUCTURES_SOLVERSTEPS_BASE_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_BASE_HH

#include<memory>
#include<type_traits>
#include<vector>


template<typename V>
class TransitionSolverStepBase
{
  public:
  // Export some types from the PDELab Vector type
  using Base = TransitionSolverStepBase<V>;
  using Vector = V;
  using GridFunctionSpace = typename Vector::GridFunctionSpace;
  using GridView = typename GridFunctionSpace::Traits::GridViewType;
  using Grid = typename GridView::Traits::Grid;
  using ctype = typename Grid::ctype;
  using Entity = typename GridView::template Codim<0>::Entity;
  using GlobalCoordinate = typename Entity::Geometry::GlobalCoordinate;
  using Range = typename Vector::field_type;
  using ConstraintsContainer = typename GridFunctionSpace::template ConstraintsContainer<Range>::Type;

  // The virtual interface - pretty simple
  virtual ~TransitionSolverStepBase() {}

  virtual void apply(Vector& vector, ConstraintsContainer& cc)
  {}
};


template<typename Vector, typename... Params>
class ParametrizedTransitionStepBase
  : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  virtual ~ParametrizedTransitionStepBase() {}

  virtual void update_transition_value(Params... val)
  {}
};


template<typename Vector, typename... Params>
class NoopParametrizationWrapper
  : public ParametrizedTransitionStepBase<Vector, Params...>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  NoopParametrizationWrapper(std::shared_ptr<TransitionSolverStepBase<Vector>> step)
    : step(step)
  {}

  virtual ~NoopParametrizationWrapper() {}

  virtual void apply(Vector& vector, typename Base::ConstraintsContainer& cc) override
  {
    step->apply(vector, cc);
  }

  private:
  std::shared_ptr<TransitionSolverStepBase<Vector>> step;
};


template<typename Vector, typename... Params>
class VariationTransitionStepBase
  : public ParametrizedTransitionStepBase<Vector, Params...>
{
  public:
  using Base = ParametrizedTransitionStepBase<Vector, Params...>;

  VariationTransitionStepBase()
    : steps(0)
  {}

  virtual ~VariationTransitionStepBase() {}

  virtual void update_transition_value(Params... val) override
  {
    for (auto step : steps)
      step->update_transition_value(val...);
  }

  template<typename STEP>
  void add(std::shared_ptr<STEP> step)
  {
    if constexpr (std::is_convertible<STEP*, ParametrizedTransitionStepBase<Vector, Params...>*>::value)
      steps.push_back(step);
    else
      steps.push_back(std::make_shared<NoopParametrizationWrapper<Vector, Params...>>(step));
  }

  template<typename STEP>
  void add(STEP& step)
  {
    add(Dune::stackobject_to_shared_ptr(step));
  }

  protected:
  std::vector<std::shared_ptr<ParametrizedTransitionStepBase<Vector, Params...>>> steps;
};


#endif
