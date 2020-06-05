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
#include<dune/structures/solversteps/adaptivity.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/constraints.hh>
#include<dune/structures/solversteps/control.hh>
#include<dune/structures/solversteps/elasticity.hh>
#include<dune/structures/solversteps/filelogger.hh>
#include<dune/structures/solversteps/interpolation.hh>
#include<dune/structures/solversteps/linearsolver.hh>
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
#include<tuple>
#include<vector>


template<typename Vector>
struct VectorToConstraintsContainer
{
  using type = typename Vector::GridFunctionSpace::template ConstraintsContainer<typename Vector::field_type>::Type;
};


template<typename... V>
class TransitionSolver
{
  public:
  using StepBase = TransitionSolverStepBase<V...>;
  using StepTraits = SimpleStepTraits<V...>;
  using EntitySet = typename StepTraits::EntitySet;
  using Parameter = typename StepTraits::Parameter;
  using Grid = typename StepTraits::Grid;

  TransitionSolver(std::shared_ptr<Grid> grid, EntitySet es)
    : grid(grid), es(es), steps(0)
  {}

  TransitionSolver(std::shared_ptr<Grid> grid, EntitySet es, std::vector<std::shared_ptr<StepBase>> steps)
    : grid(grid)
    , es(es)
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

  void apply()
  {
    // Make sure that all initial parameters are propagated to all steps.
    // This avoids dependencies on the step insertion order
    for (auto [name, param] : paramdata)
    {
      for (auto step : steps)
        step->update_parameter(name, param);
    }

    for (auto step : steps)
      step->pre();

    for (auto step : steps)
      step->apply();

    for (auto step : steps)
     step->post();
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

  std::shared_ptr<Grid> getGrid() const
  {
    return grid;
  }

  void setVectors(const std::tuple<std::shared_ptr<V>...>& vectors_)
  {
    vectors = vectors_;
  }

  void setConstraintsContainers(const std::tuple<std::shared_ptr<typename VectorToConstraintsContainer<V>::type>...>& containers_)
  {
    constraints_containers = containers_;
  }

  template<std::size_t i=0>
  auto getVector()
  {
    return std::get<i>(vectors);
  }

  auto getVectors()
  {
    return vectors;
  }

  template<std::size_t i=0>
  auto getConstraintsContainer()
  {
    return std::get<i>(constraints_containers);
  }

  private:
  std::shared_ptr<Grid> grid;
  EntitySet es;
  std::vector<std::shared_ptr<StepBase>> steps;
  std::map<std::string, Parameter> paramdata;
  std::tuple<std::shared_ptr<V>...> vectors;
  std::tuple<std::shared_ptr<typename VectorToConstraintsContainer<V>::type>...> constraints_containers;
};

#endif
