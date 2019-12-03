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


template<typename Vector, typename ValueType>
class ParametrizedMaterialStepBase
  : public ParametrizedTransitionStepBase<Vector, ValueType>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  ParametrizedMaterialStepBase(std::shared_ptr<MaterialDependantStepBase<Vector>> step,
                               std::function<void(Dune::ParameterTree&, ValueType)> modificator)
    : step(step), modificator(modificator)
  {}

  virtual ~ParametrizedMaterialStepBase() {}

  virtual void update_transition_value(ValueType val) override
  {
    auto& params = step->get_params();
    modificator(params, val);
    step->rebuild_material();
  }

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    step->apply(vector, cc);
  }

  protected:

  private:
  std::shared_ptr<MaterialDependantStepBase<Vector>> step;
  std::function<void(Dune::ParameterTree&, ValueType)> modificator;
};


template<typename Vector, typename ValueType>
class DiscreteMaterialVariationTransitionStep
  : public DiscreteVariationTransitionStep<Vector, ValueType>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  DiscreteMaterialVariationTransitionStep(
      std::function<void(Dune::ParameterTree&, ValueType)> modificator,
      std::vector<ValueType> values)
    : DiscreteVariationTransitionStep<Vector, ValueType>(values)
    , modificator(modificator)
  {}

  template<typename STEP>
  void add(std::shared_ptr<STEP> step)
  {
    if constexpr (std::is_convertible<STEP*, MaterialDependantStepBase<Vector>*>::value)
      this->steps.push_back(std::make_shared<ParametrizedMaterialStepBase<Vector, ValueType>>(step, modificator));
    else
      this->steps.push_back(std::make_shared<NoopParametrizationWrapper<Vector, ValueType>>(step));
  }

  template<typename STEP>
  void add(STEP& step)
  {
    add(Dune::stackobject_to_shared_ptr(step));
  }

  private:
  std::function<void(Dune::ParameterTree&, ValueType)> modificator;
};


template<typename Vector>
class ContinuousMaterialVariationTransitionStep
  : public ContinuousVariationTransitionStep<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  ContinuousMaterialVariationTransitionStep(std::function<void(Dune::ParameterTree&, double)> modificator, int iterations=5, double start=0.0, double end=1.0)
    : ContinuousVariationTransitionStep<Vector>(iterations, start, end)
    , modificator(modificator)
  {}

  template<typename STEP>
  void add(std::shared_ptr<STEP> step)
  {
    if constexpr (std::is_convertible<STEP*, MaterialDependantStepBase<Vector>*>::value)
    {
      std::cout << "I am adding a thing" << std::endl;
      this->steps.push_back(std::make_shared<ParametrizedMaterialStepBase<Vector, double>>(step, modificator));
    }
    else
      this->steps.push_back(std::make_shared<NoopParametrizationWrapper<Vector, double>>(step));
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
