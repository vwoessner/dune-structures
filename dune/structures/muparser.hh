#ifndef DUNE_STRUCTURES_MUPARSER_HH
#define DUNE_STRUCTURES_MUPARSER_HH

#include<dune/common/fvector.hh>
#include<dune/structures/transitionsolver.hh>
#include<muParser.h>
#include<string>


template<typename Signature>
class MuParserCallable
{};


template<typename R, typename D>
class MuParserCallable<R(D)>
{
  public:
  using Range = R;
  using Domain = D;

  MuParserCallable(std::string expr)
    : position(std::make_shared<D>(0.0))
  {
    parser.DefineVar("x", &(*position)[0]);
    parser.DefineVar("y", &(*position)[1]);
    parser.DefineVar("z", &(*position)[2]);

    parser.SetExpr(expr);
  }

  template<typename Vector>
  MuParserCallable(TransitionSolver<Vector>& solver, std::string expr)
    : MuParserCallable(expr)
  {
    for (auto [name, pointer] : solver.export_parameters())
      parser.DefineVar(name, pointer);
  }

  R operator()(const D& x)
  {
    *position = x;
    return parser.Eval();
  }

  protected:
  mu::Parser parser;
  std::shared_ptr<D> position;
};


template<typename Signature>
class MuParserTransformation
{};


template<typename R>
class MuParserTransformation<R(R, R)>
  : public MuParserCallable<R(R)>
{
  public:
  using Range = R;
  using Domain = R;

  MuParserTransformation(std::string expr)
    : MuParserCallable<R(R)>(expr)
    , solution(std::make_shared<Domain>(0.0))
  {
    define_solution_variables();
  }

  template<typename Vector>
  MuParserTransformation(TransitionSolver<Vector>& solver, std::string expr)
    : MuParserCallable<R(R)>(solver, expr)
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

#endif
