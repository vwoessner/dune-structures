#ifndef DUNE_STRUCTURES_SOLVERSTEPS_MATERIAL_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_MATERIAL_HH

#include<dune/structures/material.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/newton.hh>
#include<dune/structures/solversteps/variation.hh>

#include<functional>
#include<memory>


template<typename Vector>
class MaterialDependantStepBase
  : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  MaterialDependantStepBase(typename Base::EntitySet es,
                            std::shared_ptr<std::vector<int>> physical,
                            const Dune::ParameterTree& params
                            )
    : es(es)
    , physical(physical)
    , params(params)
    , material(parse_material<double>(es, physical, params))
  {}

  virtual ~MaterialDependantStepBase() {}

  void rebuild_material()
  {
    *material = *parse_material<double>(es, physical, params);
  }

  Dune::ParameterTree& get_params()
  {
    return params;
  }

  private:
  typename Base::EntitySet es;
  std::shared_ptr<std::vector<int>> physical;
  Dune::ParameterTree params;

  protected:
  std::shared_ptr<MaterialCollection<typename Base::EntitySet, double>> material;
};


template<typename Vector>
class ParametrizedMaterialStepBase
  : public WrapperStep<Vector, MaterialDependantStepBase<Vector>>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  ParametrizedMaterialStepBase(std::shared_ptr<MaterialDependantStepBase<Vector>> step,
                               std::function<void(Dune::ParameterTree&, typename Base::Parameter)> modificator)
    : WrapperStep<Vector, MaterialDependantStepBase<Vector>>(step)
    , modificator(modificator)
  {}

  virtual ~ParametrizedMaterialStepBase() {}

  virtual void update_parameter(std::string name, typename Base::Parameter val) override
  {
    if (name == "material")
    {
      auto& params = this->step->get_params();
      modificator(params, val);
      this->step->rebuild_material();
    }

    this->step->update_parameter(name, val);
  }

  private:
  std::function<void(Dune::ParameterTree&, typename Base::Parameter)> modificator;
};


template<typename Vector>
class DiscreteMaterialVariationTransitionStep
  : public DiscreteVariationTransitionStep<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  DiscreteMaterialVariationTransitionStep(
      std::function<void(Dune::ParameterTree&, typename Base::Parameter)> modificator,
      std::vector<typename Base::Parameter> values)
    : DiscreteVariationTransitionStep<Vector>("material", values)
    , modificator(modificator)
  {}

  template<typename STEP>
  void add(std::shared_ptr<STEP> step)
  {
    if constexpr (std::is_convertible<STEP*, MaterialDependantStepBase<Vector>*>::value)
      this->steps.push_back(std::make_shared<ParametrizedMaterialStepBase<Vector>>(step, modificator));
    else
      this->steps.push_back(step);
  }

  template<typename STEP>
  void add(STEP& step)
  {
    add(Dune::stackobject_to_shared_ptr(step));
  }

  private:
  std::function<void(Dune::ParameterTree&, typename Base::Parameter)> modificator;
};


template<typename Vector>
class ContinuousMaterialVariationTransitionStep
  : public ContinuousVariationTransitionStep<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  ContinuousMaterialVariationTransitionStep(
      std::function<void(Dune::ParameterTree&, typename Base::Parameter)> modificator,
      int iterations=5,
      double start=0.0,
      double end=1.0)
    : ContinuousVariationTransitionStep<Vector>("material", iterations, start, end)
    , modificator(modificator)
  {}

  template<typename STEP>
  void add(std::shared_ptr<STEP> step)
  {
    if constexpr (std::is_convertible<STEP*, MaterialDependantStepBase<Vector>*>::value)
    { this->steps.push_back(std::make_shared<ParametrizedMaterialStepBase<Vector>>(step, modificator));
    }
    else
      this->steps.push_back(step);
  }

  template<typename STEP>
  void add(STEP& step)
  {
    add(Dune::stackobject_to_shared_ptr(step));
  }

  private:
  std::function<void(Dune::ParameterTree&, double)> modificator;
};

#endif
