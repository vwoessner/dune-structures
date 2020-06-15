#ifndef DUNE_BLOCKLAB_CONSTRUCTION_MUPARSER_HH
#define DUNE_BLOCKLAB_CONSTRUCTION_MUPARSER_HH

/** An implementation of a MuParser-based callable that can be used
 *  to parse arithmetic expressions from configuration files.
 */

#include<dune/blocklab/utilities/stringsplit.hh>
#include<dune/common/hybridutilities.hh>
#include<dune/grid/concepts/intersection.hh>

#include<muParser.h>

#include<functional>
#include<string>


namespace Dune::BlockLab {

  // Forward declaration of the BlockSolver class
  template<typename P, typename V>
  class BlockSolver;

  template<typename Signature, typename... Callbacks>
  class MuParserCallable
  {};


  template<typename R, typename E, typename D>
  class MuParserCallable<R(E, D)>
  {
    public:
    using Range = R;
    using Domain = D;
    using Global = typename E::Geometry::GlobalCoordinate;

    template<typename... Callbacks>
    MuParserCallable(std::string expr, Callbacks&&... callbacks)
      : position(std::make_shared<Global>(0.0))
    {
      parser.DefineVar("x", &(*position)[0]);

      if constexpr (Global::dimension > 1)
        parser.DefineVar("y", &(*position)[1]);

      if constexpr (Global::dimension > 2)
        parser.DefineVar("z", &(*position)[2]);

      parser.SetExpr(expr);

      // Add the callback functions
      (this->parser.DefineFun(callbacks.first, callbacks.second, false), ...);
    }

    R operator()(const E& e, const D& x)
    {
      // Evaluate quantities accessible from muparser expressions
      *position = e.geometry().global(x);
      return static_cast<R>(parser.Eval());
    }

    protected:
    mu::Parser parser;
    std::shared_ptr<Global> position;
  };

  template<typename Signature, typename... Callbacks>
  std::function<Signature> muparser_callable(std::string expr, Callbacks&&... callbacks)
  {
    return MuParserCallable<Signature, Callbacks...>(expr, std::forward<Callbacks>(callbacks)...);
  }

  namespace impl {

    template<typename P, typename V>
    struct BlockSolverCallbacks
    {
      // The callback functions to be registered with muparser
      // Note that the signature needs to match muparser prototype exactly
      static double param(const char* name)
      {
	return solver->template param<double>(name);
      }

      // A pointer to a solver instance. This currently prohibits two solvers
      // of exact same type to be instantiated within one program
      static const BlockSolver<P, V>* solver;
    };

    template<typename P, typename V>
    const BlockSolver<P, V>* BlockSolverCallbacks<P, V>::solver = 0;

  }

  template<typename Signature, typename Solver>
  std::function<Signature> muparser_callable(std::string expr, std::shared_ptr<Solver> solver)
  {
    using CB = impl::BlockSolverCallbacks<typename Solver::ParameterTuple, typename Solver::VectorTuple>;
    CB::solver = solver.get();
    return MuParserCallable<Signature>(expr, std::make_pair("param", CB::param));
  }



  template<std::size_t size, typename Signature, typename... Callbacks>
  std::array<std::function<Signature>, size>
  muparser_callable_array(std::string expr, Callbacks&&... callbacks)
  {
    std::array<std::function<Signature>, size> result;

    auto exprs = string_split(expr);
    if (exprs.size() == 1)
      result.fill(muparser_callable<Signature, Callbacks...>(exprs[0], std::forward<Callbacks>(callbacks)...));
    else
      std::transform(exprs.begin(), exprs.end(), result.begin(),
		     [callbacks...](auto it){
                       return muparser_callable<Signature, Callbacks...>(it, std::forward<Callbacks>(callbacks)...);
                     });

    return std::move(result);
  }

  template<std::size_t size, typename Signature, typename Solver>
  std::array<std::function<Signature>, size>
  muparser_callable_array(std::string expr, std::shared_ptr<Solver> solver)
  {
    std::array<std::function<Signature>, size> result;

    auto exprs = string_split(expr);
    if (exprs.size() == 1)
      result.fill(muparser_callable<Signature, Solver>(exprs[0], solver));
    else
      std::transform(exprs.begin(), exprs.end(), result.begin(),
		     [solver](auto it){
                       return muparser_callable<Signature, Solver>(it, solver);
                     });

    return std::move(result);
  }

} // namespace Dune::BlockLab

#endif
