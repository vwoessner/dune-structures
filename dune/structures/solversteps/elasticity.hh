#ifndef DUNE_STRUCTURES_SOLVERSTEPS_ELASTICITY_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_ELASTICITY_HH

#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/instationary.hh>
#include<dune/structures/solversteps/material.hh>
#include<dune/structures/solversteps/newton.hh>

#include"operators/elasticity_operator.hh"
#include"operators/elastodynamics_spatial_operator.hh"
#include"operators/elastodynamics_temporal_operator.hh"


template<typename Vector>
class ElasticitySolverStep
  : public WrapperStep<Vector,
                       NewtonSolverTransitionStep<Vector, ElasticityOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                                             typename TransitionSolverStepBase<Vector>::GridFunctionSpace>>>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;
  using LocalOperator = ElasticityOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                           typename TransitionSolverStepBase<Vector>::GridFunctionSpace>;

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
    this->step->set_localoperator(std::make_shared<LocalOperator>(*gfs, *gfs, params, material));
    this->step->pre(vector, cc);
  }

  private:
  std::shared_ptr<typename Base::Material> material;
  Dune::ParameterTree params;
};


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

#endif
