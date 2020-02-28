#ifndef DUNE_STRUCTURES_SOLVERSTEPS_ELASTICITY_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_ELASTICITY_HH

#include<dune/structures/muparser.hh>
#include<dune/structures/elasticity.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/instationary.hh>
#include<dune/structures/solversteps/material.hh>
#include<dune/structures/solversteps/newton.hh>
#include<dune/structures/solversteps/traits.hh>


template<typename... V>
class ElasticitySolverStep
  : public WrapperStep<NewtonSolverTransitionStep<typename OperatorSwitch<typename VectorStepTraits<0, V...>::GridFunctionSpace,
                                                                          VectorStepTraits<0, V...>::dim>::Elasticity, V...>, V...>
{
  public:
  using Traits = VectorStepTraits<0, V...>;

  static constexpr int dim = Traits::dim;
  using LocalOperator = typename OperatorSwitch<typename Traits::GridFunctionSpace, dim>::Elasticity;

  ElasticitySolverStep(const Dune::ParameterTree& params)
    : WrapperStep<NewtonSolverTransitionStep<LocalOperator, V...>, V...>(std::make_shared<NewtonSolverTransitionStep<LocalOperator, V...>>(params))
    , params(params)
  {}

  virtual ~ElasticitySolverStep() {}

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
  {
    if (name == "material")
    {
      material = std::get<std::shared_ptr<typename Traits::Material>>(param);

      // I strongly dislike this lifetime dependencies, but Update parameter is available
      // very early right now - the LocalOperator might not even be constructed.
      auto lop = this->step->get_localoperator();
      if (lop != nullptr)
        lop->setMaterial(material);
    }

    this->step->update_parameter(name, param);
  }

  virtual void pre() override
  {
    auto vector = this->solver->getVector();
    auto gfs = vector->gridFunctionSpaceStorage();
    auto lop = std::make_shared<LocalOperator>(*gfs, *gfs, params, material);
    force = std::make_shared<typename Traits::Vector>(*gfs, 0.0);
    traction = std::make_shared<typename Traits::Vector>(*gfs, 0.0);
    lop->setCoefficientForce(gfs, force);
    lop->setCoefficientTraction(gfs, traction);
    this->step->set_localoperator(lop);
    this->step->pre();
  }

  virtual void apply() override
  {
    auto vector = this->solver->getVector();
    auto constraintscontainer = this->solver->getConstraintsContainer();
    auto gfs = vector->gridFunctionSpaceStorage();

    // Interpolate the force vector
    auto force_gf = makeGridFunctionTreeFromCallables(*gfs, get_callable_array<double(typename Traits::Entity, typename Traits::GlobalCoordinate), V...>(*(this->solver), params.get<std::string>("force", "0.0")));
    Dune::PDELab::interpolate(force_gf, *gfs, *force);

    // Interpolate the traction vector
    auto traction_gf = makeGridFunctionTreeFromCallables(*gfs, get_callable_array<double(typename Traits::Entity, typename Traits::GlobalCoordinate), V...>(*(this->solver), params.get<std::string>("traction", "0.0")));
    Dune::PDELab::interpolate(traction_gf, *gfs, *traction);

    this->step->apply();
  }

  private:
  std::shared_ptr<typename Traits::Material> material;
  std::shared_ptr<typename Traits::Vector> force;
  std::shared_ptr<typename Traits::Vector> traction;
  Dune::ParameterTree params;
};


template<typename... V>
class QuasiStaticElastoDynamicsSolverStep
  : public WrapperStep<OneStepMethodStep<typename OperatorSwitch<typename VectorStepTraits<0, V...>::GridFunctionSpace,
                                                                 VectorStepTraits<0, V...>::dim>::Elasticity,
                                         typename OperatorSwitch<typename VectorStepTraits<0, V...>::GridFunctionSpace,
                                                                 VectorStepTraits<0, V...>::dim>::Mass,
                                         V...>, V...>
{
  public:
  using Traits = VectorStepTraits<0, V...>;

  using SpatialLocalOperator = typename OperatorSwitch<typename Traits::GridFunctionSpace,
                                                       Traits::dim>::Elasticity;

  using TemporalLocalOperator = typename OperatorSwitch<typename Traits::GridFunctionSpace,
                                                        Traits::dim>::Mass;

  QuasiStaticElastoDynamicsSolverStep(const Dune::ParameterTree& rootparams)
    : WrapperStep<OneStepMethodStep<SpatialLocalOperator, TemporalLocalOperator, V...>, V...>(std::make_shared<OneStepMethodStep<SpatialLocalOperator, TemporalLocalOperator, V...>>())
    , params(rootparams)
  {}

  virtual ~QuasiStaticElastoDynamicsSolverStep() override {}

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
  {
    if (name == "material")
    {
      material = std::get<std::shared_ptr<typename Traits::Material>>(param);

      // I strongly dislike this lifetime dependencies, but Update parameter is available
      // very early right now - the LocalOperator might not even be constructed.
      auto lop = this->step->get_spatial_localoperator();
      if (lop != nullptr)
        lop->setMaterial(material);
    }

    this->step->update_parameter(name, param);
  }

  virtual void pre() override
  {
    auto vector = this->solver->getVector();
    auto& gfs = vector->gridFunctionSpace();
    this->step->set_spatial_localoperator(std::make_shared<SpatialLocalOperator>(gfs, gfs, params, material));
    this->step->set_temporal_localoperator(std::make_shared<TemporalLocalOperator>(gfs, gfs, params));
    this->step->pre();
  }

  private:
  std::shared_ptr<typename Traits::Material> material;
  Dune::ParameterTree params;
};


#endif
