#ifndef DUNE_STRUCTURES_SOLVERSTEPS_MATERIAL_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_MATERIAL_HH

#include<dune/structures/material.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/traits.hh>
#include<dune/structures/solversteps/variation.hh>

#include<functional>
#include<memory>


template<typename... V>
class MaterialInitialization
  : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  MaterialInitialization(typename Traits::EntitySet es,
      std::shared_ptr<std::vector<int>> physical,
      const Dune::ParameterTree& params,
      const Dune::ParameterTree& rootparams
      )
    : es(es)
    , physical(physical)
    , params(params)
    , rootparams(rootparams)
    , material(parse_material<double>(es, physical, params, rootparams))
  {}

  virtual ~MaterialInitialization() {}

  virtual void set_solver(std::shared_ptr<typename Traits::Solver> solver_) override
  {
    this->solver = solver_;
    this->solver->introduce_parameter("physical", physical);
    this->solver->introduce_parameter("material_params", params);
    this->solver->introduce_parameter("material", material);
  }

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
  {
    if (name == "material_params")
    {
      params = std::get<Dune::ParameterTree>(param);
      this->solver->update_parameter("material", parse_material<double>(es, physical, params, rootparams));
    }
  }

  protected:
  typename Traits::EntitySet es;
  std::shared_ptr<std::vector<int>> physical;
  Dune::ParameterTree params;
  Dune::ParameterTree rootparams;
  std::shared_ptr<typename Traits::Material> material;
  std::shared_ptr<typename Traits::Solver> solver;
};


template<typename... V>
class DiscreteMaterialVariationTransitionStep
  : public DiscreteVariationTransitionStep<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  DiscreteMaterialVariationTransitionStep(
    std::string pname,
    std::function<void(Dune::ParameterTree&, typename Traits::Parameter)> modificator,
    std::vector<typename Traits::Parameter> values)
    : DiscreteVariationTransitionStep<V...>(pname, values)
    , pname(pname)
    , modificator(modificator)
  {}

  virtual ~DiscreteMaterialVariationTransitionStep() {}

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
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
  std::function<void(Dune::ParameterTree&, typename Traits::Parameter)> modificator;
};


template<typename... V>
class ContinuousMaterialVariationTransitionStep
  : public ContinuousVariationTransitionStep<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  ContinuousMaterialVariationTransitionStep(
    std::string pname,
    std::function<void(Dune::ParameterTree&, typename Traits::Parameter)> modificator,
    int iterations=5,
    double start=0.0,
    double end=1.0)
    : ContinuousVariationTransitionStep<V...>(pname, iterations, start, end)
    , pname(pname)
    , modificator(modificator)
  {}

  virtual ~ContinuousMaterialVariationTransitionStep() {}

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
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
