#ifndef DUNE_STRUCTURES_SOLVERSTEPS_ELASTICITY_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_ELASTICITY_HH

#include<dune/structures/muparser.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/instationary.hh>
#include<dune/structures/solversteps/material.hh>
#include<dune/structures/solversteps/newton.hh>

#include"operators/elasticity_operator.hh"
#include"operators/quasistatic_mass_operator.hh"
//#include"operators/elastodynamics_spatial_operator.hh"
//#include"operators/elastodynamics_temporal_operator.hh"


template<typename Vector>
class ElasticitySolverStep
  : public WrapperStep<Vector,
                       NewtonSolverTransitionStep<Vector, ElasticityOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                                             typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                                             typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                                             typename TransitionSolverStepBase<Vector>::GridFunctionSpace>>>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;
  using LocalOperator = ElasticityOperator<typename Base::GridFunctionSpace,
                                           typename Base::GridFunctionSpace,
                                           typename Base::GridFunctionSpace,
                                           typename Base::GridFunctionSpace>;

  ElasticitySolverStep(const Dune::ParameterTree& params)
    : WrapperStep<Vector, NewtonSolverTransitionStep<Vector, LocalOperator>>(std::make_shared<NewtonSolverTransitionStep<Vector, LocalOperator>>())
    , params(params)
  {}

  virtual ~ElasticitySolverStep() {}

  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    if (name == "material")
    {
      material = std::get<std::shared_ptr<typename Base::Material>>(param);

      // I strongly dislike this lifetime dependencies, but Update parameter is available
      // very early right now - the LocalOperator might not even be constructed.
      auto lop = this->step->get_localoperator();
      if (lop != nullptr)
        lop->setMaterial(material);
    }

    this->step->update_parameter(name, param);
  }

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    auto gfs = vector->gridFunctionSpaceStorage();
    auto lop = std::make_shared<LocalOperator>(*gfs, *gfs, params, material);
    force = std::make_shared<Vector>(*gfs, 0.0);
    traction = std::make_shared<Vector>(*gfs, 0.0);
    lop->setCoefficientForce(gfs, force);
    lop->setCoefficientTraction(gfs, traction);
    this->step->set_localoperator(lop);
    this->step->pre(vector, cc);
  }

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    auto gfs = vector->gridFunctionSpaceStorage();

    // Interpolate the force vector
    auto force_gf = makeGridFunctionTreeFromCallables(*gfs, get_callable_array<Vector, double(typename Base::Entity, typename Base::GlobalCoordinate)>(*(this->solver), params.get<std::string>("force", "0.0")));
    Dune::PDELab::interpolate(force_gf, *gfs, *force);

    // Interpolate the traction vector
    auto traction_gf = makeGridFunctionTreeFromCallables(*gfs, get_callable_array<Vector, double(typename Base::Entity, typename Base::GlobalCoordinate)>(*(this->solver), params.get<std::string>("traction", "0.0")));
    Dune::PDELab::interpolate(traction_gf, *gfs, *traction);

    this->step->apply(vector, cc);
  }

  private:
  std::shared_ptr<typename Base::Material> material;
  std::shared_ptr<Vector> force;
  std::shared_ptr<Vector> traction;
  Dune::ParameterTree params;
};


template<typename Vector>
class QuasiStaticElastoDynamicsSolverStep
  : public WrapperStep<Vector, OneStepMethodStep<Vector,
                                                 ElasticityOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                                    typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                                    typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                                    typename TransitionSolverStepBase<Vector>::GridFunctionSpace>,
                                                 QuasiStaticMassOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                                         typename TransitionSolverStepBase<Vector>::GridFunctionSpace>
    >>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  using SpatialLocalOperator = ElasticityOperator<typename Base::GridFunctionSpace,
                                                  typename Base::GridFunctionSpace,
                                                  typename Base::GridFunctionSpace,
                                                  typename Base::GridFunctionSpace>;

  using TemporalLocalOperator = QuasiStaticMassOperator<typename Base::GridFunctionSpace,
                                                        typename Base::GridFunctionSpace>;

  QuasiStaticElastoDynamicsSolverStep(const Dune::ParameterTree& rootparams)
    : WrapperStep<Vector, OneStepMethodStep<Vector, SpatialLocalOperator, TemporalLocalOperator>>(std::make_shared<OneStepMethodStep<Vector, SpatialLocalOperator, TemporalLocalOperator>>())
    , params(rootparams)
  {}

  virtual ~QuasiStaticElastoDynamicsSolverStep() override {}

  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    if (name == "material")
    {
      material = std::get<std::shared_ptr<typename Base::Material>>(param);

      // I strongly dislike this lifetime dependencies, but Update parameter is available
      // very early right now - the LocalOperator might not even be constructed.
      auto lop = this->step->get_spatial_localoperator();
      if (lop != nullptr)
        lop->setMaterial(material);
    }

    this->step->update_parameter(name, param);
  }

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    auto& gfs = vector->gridFunctionSpace();
    this->step->set_spatial_localoperator(std::make_shared<SpatialLocalOperator>(gfs, gfs, params, material));
    this->step->set_temporal_localoperator(std::make_shared<TemporalLocalOperator>(gfs, gfs, params));
    this->step->pre(vector, cc);
  }

  private:
  std::shared_ptr<typename Base::Material> material;
  Dune::ParameterTree params;
};


/*
template<typename Vector>
class ElastoDynamicsSolverStep
  :  public WrapperStep<Vector, OneStepMethodStep<Vector,
                                         ElastoDynamicsSpatialOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                                       typename TransitionSolverStepBase<Vector>::GridFunctionSpace>,
                                         ElastoDynamicsTemporalOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                                        typename TransitionSolverStepBase<Vector>::GridFunctionSpace>
                                         >>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  using SpatialLocalOperator = ElastoDynamicsSpatialOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                             typename TransitionSolverStepBase<Vector>::GridFunctionSpace>;

  using TemporalLocalOperator = ElastoDynamicsTemporalOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                               typename TransitionSolverStepBase<Vector>::GridFunctionSpace>;

  ElastoDynamicsSolverStep(const Dune::ParameterTree& rootparams)
    : WrapperStep<Vector, OneStepMethodStep<Vector, SpatialLocalOperator, TemporalLocalOperator>>(std::make_shared<OneStepMethodStep<Vector, SpatialLocalOperator, TemporalLocalOperator>>())
    , params(rootparams)
  {}

  template<typename... FUNCS>
  ElastoDynamicsSolverStep(const Dune::ParameterTree& rootparams,
                           const FUNCS&... funcs)
   : WrapperStep<Vector, OneStepMethodStep<Vector, SpatialLocalOperator, TemporalLocalOperator>>(std::make_shared<VariableBoundaryOneStepMethodStep<Vector, SpatialLocalOperator, TemporalLocalOperator>>(funcs...))
  , params(rootparams)
  {}

  template<typename... FUNCS>
  ElastoDynamicsSolverStep(const Dune::ParameterTree& rootparams,
                           const std::array<std::function<typename VariableBoundaryOneStepMethodStep<Vector, SpatialLocalOperator, TemporalLocalOperator>::FunctionSignature>,
                                            Dune::TypeTree::TreeInfo<typename Base::GridFunctionSpace>::leafCount>& funcs)
   : WrapperStep<Vector, OneStepMethodStep<Vector, SpatialLocalOperator, TemporalLocalOperator>>(std::make_shared<VariableBoundaryOneStepMethodStep<Vector, SpatialLocalOperator, TemporalLocalOperator>>(funcs))
   , params(rootparams)
  {}

  virtual ~ElastoDynamicsSolverStep() {}

  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    if (name == "material")
    {
      material = std::get<std::shared_ptr<typename Base::Material>>(param);

      // I strongly dislike this lifetime dependencies, but Update parameter is available
      // very early right now - the LocalOperator might not even be constructed.
      auto lop = this->step->get_spatial_localoperator();
      if (lop != nullptr)
        lop->setMaterial(material);
    }

    this->step->update_parameter(name, param);
  }

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    auto& gfs = vector->gridFunctionSpace();
    this->step->set_spatial_localoperator(std::make_shared<SpatialLocalOperator>(gfs, gfs, params, material));
    this->step->set_temporal_localoperator(std::make_shared<TemporalLocalOperator>(gfs, gfs, params));
    this->step->pre(vector, cc);
  }

  private:
  std::shared_ptr<typename Base::Material> material;
  Dune::ParameterTree params;
};
*/

#endif
