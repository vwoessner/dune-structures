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
  {
    parser.DefineVar("x", &position[0]);
    parser.DefineVar("y", &position[1]);
    parser.DefineVar("z", &position[2]);

    parser.SetExpr(expr);
  }

  template<typename Vector>
  MuParserCallable(const TransitionSolver<Vector>& solver, std::string expr)
    : MuParserCallable(expr)
  {
    for (auto [name, pointer] : solver.export_parameters())
      parser.DefineVar(name, pointer);
  }

  R operator()(const D& x)
  {
    position = x;
    return parser.Eval();
  }

  private:
  mu::Parser parser;
  D position;
};


#endif
