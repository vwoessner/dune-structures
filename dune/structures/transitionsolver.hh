#ifndef DUNE_STRUCTURES_TRANSITIONSOLVER_HH
#define DUNE_STRUCTURES_TRANSITIONSOLVER_HH

/**
 *  A nonlinear solver that solves for a sequence of problems
 *  of increasing physical complexity
 */
#include<dune/common/shared_ptr.hh>
#include<dune/pdelab.hh>
#include<dune/structures/utilities.hh>

// Include all available solver steps in order to have this header serve
// as a convenience header
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/constraints.hh>
#include<dune/structures/solversteps/elasticity.hh>
#include<dune/structures/solversteps/interpolation.hh>
#include<dune/structures/solversteps/material.hh>
#include<dune/structures/solversteps/newton.hh>
#include<dune/structures/solversteps/onetoone.hh>
#include<dune/structures/solversteps/parameter.hh>
#include<dune/structures/solversteps/probe.hh>
#include<dune/structures/solversteps/reinforced.hh>
#include<dune/structures/solversteps/transformation.hh>
#include<dune/structures/solversteps/variation.hh>
#include<dune/structures/solversteps/visualization.hh>

#include<map>
#include<memory>
#include<vector>


template<typename... V>
class TransitionSolver
{
  public:
  using StepBase = TransitionSolverStepBase<V...>;
  using StepTraits = SimpleStepTraits<V...>;
  using EntitySet = typename StepTraits::EntitySet;
  using ConstraintsContainer = typename StepTraits::ConstraintsContainer;
  using Parameter = typename StepTraits::Parameter;

  TransitionSolver(EntitySet es)
    : es(es), steps(0)
  {}

  TransitionSolver(EntitySet es, std::vector<std::shared_ptr<StepBase>> steps)
    : es(es)
    , steps(steps)
  {}

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

  void apply(std::shared_ptr<typename StepTraits::Vector> vector, std::shared_ptr<ConstraintsContainer> constraintscontainer)
  {
    // Make sure that all initial parameters are propagated to all steps.
    // This avoids dependencies on the step insertion order
    for (auto [name, param] : paramdata)
    {
      for (auto step : steps)
        step->update_parameter(name, param);
    }

    for (auto step : steps)
      step->pre(vector, constraintscontainer);

    for (auto step : steps)
      step->apply(vector, constraintscontainer);

    for (auto step : steps)
     step->post(vector, constraintscontainer);
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
    std::get<typename std::decay<T>::type>(paramdata[name]) = std::forward<T>(val);

    for (auto step: steps)
      step->update_parameter(name, paramdata[name]);
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

  // Additional parameter interface for the time
  // PDELab expects this type of interface from a time-tracking object
  double getTime() const
  {
    return param<double>("time");
  }

  // Additional parameter interface for the time
  // PDELab expects this type of interface from a time-tracking object
  void setTime(double t)
  {
    update_parameter("time", t);
  }

  std::vector<std::pair<std::string, double*>> export_parameters()
  {
    std::vector<std::pair<std::string, double*>> ret;
    for (auto it = paramdata.begin(); it != paramdata.end(); ++it)
    {
      double* val = std::get_if<double>(&it->second);
      if (val != 0)
        ret.push_back(std::make_pair(it->first, val));
    }

    return ret;
  }

  bool has_parameter(std::string name)
  {
    return paramdata.count(name) > 0;
  }

  EntitySet entitySet() const
  {
    return es;
  }

  private:
  EntitySet es;
  std::vector<std::shared_ptr<StepBase>> steps;
  std::map<std::string, Parameter> paramdata;
};

#endif
