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

  // The possible types for parametrization of solver steps
  using Material = MaterialCollection<EntitySet, double>;
  using Parameter = std::variant<bool,
                                 double,
                                 int,
                                 std::string,
                                 Dune::FieldVector<double, dim>,
                                 std::shared_ptr<Material>,
                                 std::shared_ptr<std::vector<int>>,
                                 Dune::ParameterTree,
                                 std::shared_ptr<AbstractLocalOperatorInterface<typename V::GridFunctionSpace>>...>;

};


template<std::size_t i, typename... V>
class VectorStepTraits
  : public SimpleStepTraits<V...>
{
  using Base = SimpleStepTraits<V...>;

  public:
  // Export all types from the simple traits class
  using Base::Solver;
  using Base::GridView;
  using Base::Grid;
  using Base::EntitySet;
  using Base::ctype;
  using Base::Entity;
  using Base::GlobalCoordinate;
  using Base::Range;
  using Base::Material;
  using Base::Parameter;

  // Export the vector dependent types
  using Vector = typename std::tuple_element<i, std::tuple<V...>>::type;
  using GridFunctionSpace = typename Vector::GridFunctionSpace;
  using ConstraintsContainer = typename GridFunctionSpace::template ConstraintsContainer<typename Base::Range>::Type;
  using VectorBackend = typename GridFunctionSpace::Traits::Backend;
};

#endif
