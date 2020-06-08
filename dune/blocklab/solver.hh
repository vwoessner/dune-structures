#ifndef DUNE_BLOCKLAB_SOLVER_HH
#define DUNE_BLOCKLAB_SOLVER_HH

/** This is the main include for the dune-blocklab solver class.
 *  No additional includes needed.
 */

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/blocks/blocktraits.hh>
#include<dune/common/shared_ptr.hh>

#include<memory>
#include<string>
#include<vector>

namespace Dune::BlockLab {

  template<typename P, typename V>
  class BlockSolver
  {
    public:
    using Traits = BlockTraits<P, V>;
    using Parameter = typename Traits::Parameter;

    void apply()
    {
      for (auto step : steps)
	step->setup();

      for (auto step : steps)
	step->apply();
    }

    template<typename STEP>
    void add(std::shared_ptr<STEP> step)
    {
      steps.push_back(step);
      step->set_solver(Dune::stackobject_to_shared_ptr(*this));
    }

    template<typename STEP>
    void add(STEP& step)
    {
      add(Dune::stackobject_to_shared_ptr(step));
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

    void update_parameter(std::string name, Parameter& val)
    {
      paramdata[name] = val;

      for (auto step: steps)
        step->update_parameter(name, paramdata[name]);
    }

    template<typename T>
    T param(std::string name) const
    {
      return std::get<typename std::decay<T>::type>(paramdata.find(name)->second);
    }

    private:
    std::vector<std::shared_ptr<AbstractBlockBase<P, V>>> steps;
    std::map<std::string, Parameter> paramdata;
  };

} // namespace Dune::BlockLab

#endif
