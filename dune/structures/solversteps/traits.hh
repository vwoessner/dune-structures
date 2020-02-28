#ifndef DUNE_STRUCTURES_SOLVERSTEPS_TRAITS_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_TRAITS_HH

#include<dune/common/parametertree.hh>
#include<dune/structures/material.hh>

#include<memory>
#include<string>
#include<tuple>
#include<variant>
#include<vector>


// Forward declaration of the solver class
template<typename... V>
class TransitionSolver;


template<typename... V>
class SimpleStepTraits
{
  // Extracting the first vector and use its grid information
  using First = typename std::tuple_element<0, std::tuple<V...>>::type;

  public:
  // The transition solver for this
  using Solver = TransitionSolver<V...>;

  // Export some grid related types
  using GridView = typename First::GridFunctionSpace::Traits::GridViewType;
  static constexpr int dim = GridView::dimension;
  using Grid = typename GridView::Traits::Grid;
  using EntitySet = typename First::GridFunctionSpace::Traits::EntitySet;
  using ctype = typename Grid::ctype;
  using Entity = typename GridView::template Codim<0>::Entity;
  using GlobalCoordinate = typename Entity::Geometry::GlobalCoordinate;
  using Range = typename First::field_type;

  // These need to go away in the transition
  using Vector = First;
  using GridFunctionSpace = typename Vector::GridFunctionSpace;
  using ConstraintsContainer = typename GridFunctionSpace::template ConstraintsContainer<Range>::Type;
  using VectorBackend = typename GridFunctionSpace::Traits::Backend;

  // The possible types for parametrization of solver steps
  using Material = MaterialCollection<EntitySet, double>;
  using Parameter = std::variant<bool,
                                 double,
                                 int,
                                 std::string,
                                 std::shared_ptr<Material>,
                                 std::shared_ptr<std::vector<int>>,
                                 Dune::ParameterTree>;

};

#endif
