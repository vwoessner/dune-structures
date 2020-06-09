#ifndef DUNE_BLOCKLAB_SOLVER_HH
#define DUNE_BLOCKLAB_SOLVER_HH

/** This is the main include for the dune-blocklab solver class.
 *  No additional includes needed.
 */

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/blocks/blocktraits.hh>
#include<dune/blocklab/operators/virtualinterface.hh>
#include<dune/blocklab/utilities/uniquevariant.hh>
#include<dune/common/fvector.hh>
#include<dune/common/shared_ptr.hh>

#include<map>
#include<memory>
#include<string>
#include<tuple>
#include<vector>

namespace Dune::BlockLab {

  template<typename P, typename V>
  class BlockSolver;

  template<typename... P, typename... V>
  class BlockSolver<std::tuple<P...>, std::tuple<V...>>
  {
    public:
    // Export the template parameters of the solver class
    using ParameterTuple = std::tuple<P...>;
    using VectorTuple = std::tuple<V...>;

    // The possible types for parametrization of solver steps
    using Parameter = unique_variant<bool,
                                     double,
                                     int,
                                     std::string,
                                     Dune::ParameterTree,
                                     std::shared_ptr<AbstractLocalOperatorInterface<typename V::GridFunctionSpace>>...,
                                     P...>;

    void apply()
    {
      for (auto block : blocks)
	block->setup();

      for (auto block : blocks)
	block->apply();
    }

    template<typename BLOCK>
    void add(std::shared_ptr<BLOCK> block)
    {
      blocks.push_back(block);
      block->set_solver(Dune::stackobject_to_shared_ptr(*this));
    }

    template<typename BLOCK>
    void add(BLOCK& block)
    {
      add(Dune::stackobject_to_shared_ptr(block));
    }

    template<typename T>
    void introduce_parameter(std::string name, T&& val)
    {
      introduce_parameter(name, Parameter(std::forward<T>(val)));
    }

    void introduce_parameter(std::string name, Parameter param)
    {
      if (paramdata.count(name) == 0)
        paramdata[name] = param;
    }

    template<typename T>
    void update_parameter(std::string name, T&& val)
    {
      update_parameter(name, Parameter(std::forward<T>(val)));
    }

    void update_parameter(std::string name, Parameter val)
    {
      paramdata[name] = val;

      for (auto block: blocks)
        block->update_parameter(name, paramdata[name]);
    }

    template<typename T>
    T param(std::string name) const
    {
      return std::get<typename std::decay<T>::type>(paramdata.find(name)->second);
    }

    private:
    std::vector<std::shared_ptr<AbstractBlockBase<std::tuple<P...>, std::tuple<V...>>>> blocks;
    std::map<std::string, Parameter> paramdata;
  };

} // namespace Dune::BlockLab

#endif
