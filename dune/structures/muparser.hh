#ifndef DUNE_STRUCTURES_MUPARSER_HH
#define DUNE_STRUCTURES_MUPARSER_HH

#include<dune/common/fvector.hh>
#include<dune/grid/concepts/intersection.hh>
#include<dune/structures/utilities.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/traits.hh>

#include<muParser.h>
#include<string>


// Forward declaration of the transition solver class
template<typename... V>
class TransitionSolver;


template<typename Signature, typename... Vector>
class MuParserCallable
{};


template<typename R, typename E, typename D, typename... Vector>
class MuParserCallable<R(E, D), Vector...>
{
  public:
  using Range = R;
  using Domain = D;
  using Global = typename E::Geometry::GlobalCoordinate;

  MuParserCallable(TransitionSolver<Vector...>& solver, std::string expr)
    : position(std::make_shared<Global>(0.0))
    , current_physical(std::make_shared<double>(0))
    , indexset(solver.entitySet().indexSet())
    , physical_info(solver.template param<std::shared_ptr<std::vector<int>>>("physical"))
  {
    parser.DefineVar("x", &(*position)[0]);
    parser.DefineVar("y", &(*position)[1]);
    parser.DefineVar("z", &(*position)[2]);

    parser.SetExpr(expr);

    parser.DefineVar("group", current_physical.get());
    for (auto [name, pointer] : solver.export_parameters())
      parser.DefineVar(name, pointer);
  }

  R operator()(const E& e, const D& x)
  {
    // Look up the index of this cell (or the neighboring in the intersection case)
    int index;
    if constexpr (Dune::isIntersection<E>())
        index = indexset.index(e.inside());
    else
        index = indexset.index(e);

    // Evaluate quantities accessible from muparser expressions
    *current_physical = (*physical_info)[index];
    *position = e.geometry().global(x);
    return parser.Eval();
  }

  protected:
  mu::Parser parser;
  std::shared_ptr<Global> position;
  std::shared_ptr<double> current_physical;
  const typename TransitionSolver<Vector...>::EntitySet::Traits::IndexSet& indexset;
  std::shared_ptr<std::vector<int>> physical_info;

};


template<typename Signature, typename... Vector>
class MuParserTransformation
{};


template<typename R, typename... Vector>
class MuParserTransformation<R(R, R), Vector...>
  : public MuParserCallable<R(R), Vector...>
{
  public:
  using Range = R;
  using Domain = R;

  MuParserTransformation(TransitionSolver<Vector...>& solver, std::string expr)
    : MuParserCallable<R(R), Vector...>(solver, expr)
    , solution(std::make_shared<Domain>(0.0))
  {
    define_solution_variables();
  }

  Range operator()(const Range& sol, const Domain& x)
  {
    *this->position = x;
    *solution = sol;
    int num;
    auto eval = this->parser.Eval(num);

    return Dune::FieldVector<double, 3>{eval[0], eval[1], eval[2]};
  }

  private:
  void define_solution_variables()
  {
    this->parser.DefineVar("ux", &(*solution)[0]);
    this->parser.DefineVar("uy", &(*solution)[1]);
    this->parser.DefineVar("uz", &(*solution)[2]);
  }
  std::shared_ptr<Domain> solution;
};


template<typename Signature, typename... Vector>
std::function<Signature> get_callable(TransitionSolver<Vector...>& solver, std::string expr)
{
  return MuParserCallable<Signature, Vector...>(solver, expr);
}


template<typename Signature, typename... Vector>
std::function<Signature> get_transformation(TransitionSolver<Vector...>& solver, std::string expr)
{
  return MuParserTransformation<Signature, Vector...>(solver, expr);
}


template<typename Signature, typename... Vector>
std::array<std::function<Signature>,
           Dune::TypeTree::TreeInfo<typename VectorStepTraits<0, Vector...>::GridFunctionSpace>::leafCount> get_callable_array(TransitionSolver<Vector...>& solver, std::string expr)
{
  using GFS = typename VectorStepTraits<0, Vector...>::GridFunctionSpace;
  constexpr auto len = Dune::TypeTree::TreeInfo<GFS>::leafCount;
  std::array<std::function<Signature>, len> result;

  auto exprs = str_split(expr);
  if (exprs.size() == 1)
    result.fill(get_callable<Signature, Vector...>(solver, exprs[0]));
  else
    std::transform(exprs.begin(), exprs.end(), result.begin(), [&solver](auto it){ return get_callable<Signature, Vector...>(solver, it); });

  return std::move(result);
}

#endif
