#ifndef DUNE_STRUCTURES_SOLVERSTEPS_MATERIAL_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_MATERIAL_HH

#include<dune/structures/material.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/newton.hh>
#include<dune/structures/solversteps/variation.hh>

#include<functional>
#include<memory>


template<typename Vector>
class MaterialInitialization
  : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  MaterialInitialization(typename Base::EntitySet es,
      std::shared_ptr<std::vector<int>> physical,
      const Dune::ParameterTree& params
      )
    : es(es)
    , physical(physical)
    , params(params)
    , material(parse_material<double>(es, physical, params))
  {}

  virtual ~MaterialInitialization() {}

  virtual void set_solver(std::shared_ptr<typename Base::Solver> solver_) override
  {
    this->solver = solver_;
    this->solver->update_parameter("material_params", params);
  }

  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    if (name == "material_params")
    {
      params = std::get<Dune::ParameterTree>(param);
      this->solver->update_parameter("material", parse_material<double>(es, physical, params));
    }
  }

  protected:
  typename Base::EntitySet es;
  std::shared_ptr<std::vector<int>> physical;
  Dune::ParameterTree params;
  std::shared_ptr<typename Base::Material> material;
  std::shared_ptr<typename Base::Solver> solver;
};


template<typename Vector>
class DiscreteMaterialVariationTransitionStep
  : public DiscreteVariationTransitionStep<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  DiscreteMaterialVariationTransitionStep(
    std::string pname,
    std::function<void(Dune::ParameterTree&, typename Base::Parameter)> modificator,
    std::vector<typename Base::Parameter> values)
    : DiscreteVariationTransitionStep<Vector>(pname, values)
    , pname(pname)
    , modificator(modificator)
  {}

  virtual ~DiscreteMaterialVariationTransitionStep() {}

  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    if (pname == name)
    {
      Dune::ParameterTree config = this->solver->template param<Dune::ParameterTree>("material_params");
      modificator(config, param);
      this->solver->update_parameter("material_params", config);
    }

    for (auto step : this->steps)
      step->update_parameter(name, param);
  }

  private:
  std::string pname;
  std::function<void(Dune::ParameterTree&, typename Base::Parameter)> modificator;
};


template<typename Vector>
class ContinuousMaterialVariationTransitionStep
  : public ContinuousVariationTransitionStep<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  ContinuousMaterialVariationTransitionStep(
    std::string pname,
    std::function<void(Dune::ParameterTree&, typename Base::Parameter)> modificator,
    int iterations=5,
    double start=0.0,
    double end=1.0)
    : ContinuousVariationTransitionStep<Vector>(pname, iterations, start, end)
    , pname(pname)
    , modificator(modificator)
  {}

  virtual ~ContinuousMaterialVariationTransitionStep() {}

  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    if (pname == name)
    {
      Dune::ParameterTree config = this->solver->template param<Dune::ParameterTree>("material_params");
      modificator(config, param);
      this->solver->update_parameter("material_params", config);
    }

    for (auto step : this->steps)
      step->update_parameter(name, param);
  }

  private:
  std::string pname;
  std::function<void(Dune::ParameterTree&, double)> modificator;
};

#endif
