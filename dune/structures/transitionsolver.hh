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
#include<dune/structures/solversteps/probe.hh>
#include<dune/structures/solversteps/transformation.hh>
#include<dune/structures/solversteps/variation.hh>
#include<dune/structures/solversteps/visualization.hh>

#include<map>
#include<memory>
#include<vector>


template<typename Vector>
class TransitionSolver
{
  public:
  using ConstraintsContainer = typename Vector::GridFunctionSpace::template ConstraintsContainer<typename Vector::field_type>::Type;
  using Parameter = typename TransitionSolverStepBase<Vector>::Parameter;

  TransitionSolver()
    : steps(0)
  {}

  TransitionSolver(std::vector<std::shared_ptr<TransitionSolverStepBase<Vector>>> steps)
    : steps(steps)
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

  void apply(std::shared_ptr<Vector> vector, std::shared_ptr<ConstraintsContainer> constraintscontainer)
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

  void update_parameter(std::string name, Parameter p)
  {
    paramdata[name] = p;

    for (auto step: steps)
      step->update_parameter(name, p);
  }

  template<typename T>
  T param(std::string name) const
  {
    return std::get<T>(paramdata.find(name)->second);
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

  std::vector<std::pair<std::string, Parameter*>> export_parameters() const
  {
    std::vector<std::pair<std::string, Parameter*>> ret;
    std::transform(paramdata.begin(), paramdata.end(), ret.begin(),
                   [](auto p){ return std::make_pair<std::string, Parameter*>(p.first, &p.second); });
    return ret;
  }

  private:
  std::vector<std::shared_ptr<TransitionSolverStepBase<Vector>>> steps;
  std::map<std::string, Parameter> paramdata;
};

#endif
