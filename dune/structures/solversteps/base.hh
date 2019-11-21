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


#endif
