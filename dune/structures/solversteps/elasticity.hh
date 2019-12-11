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
  : public MaterialDependantStepBase<Vector>
  , public WrapperStep<Vector,
                       NewtonSolverTransitionStep<Vector, ElasticityOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                                             typename TransitionSolverStepBase<Vector>::GridFunctionSpace>>>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;
  using LocalOperator = ElasticityOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                           typename TransitionSolverStepBase<Vector>::GridFunctionSpace>;

  ElasticitySolverStep(typename Base::EntitySet es,
                       std::shared_ptr<std::vector<int>> physical,
                       const Dune::ParameterTree& rootparams
                       )
    : MaterialDependantStepBase<Vector>(es, physical, rootparams.sub("material"))
    , WrapperStep<Vector, NewtonSolverTransitionStep<Vector, LocalOperator>>(std::make_shared<NewtonSolverTransitionStep<Vector, LocalOperator>>())
    , params(rootparams)
  {}

  virtual ~ElasticitySolverStep() {}

  //TODO: Get rid of this which is necessary because of multiple inheritance
  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    this->step->update_parameter(name, param);
  }

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    auto gfs = vector->gridFunctionSpaceStorage();
    this->step->set_localoperator(std::make_shared<LocalOperator>(*gfs, *gfs, params, this->material));
    this->step->pre(vector, cc);
  }

  //TODO: Get rid of this which is necessary because of multiple inheritance
  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    this->step->apply(vector, cc);
  }

  private:
  Dune::ParameterTree params;
};


template<typename Vector>
class ElastoDynamicsSolverStep
  : public MaterialDependantStepBase<Vector>
  , public WrapperStep<Vector, OneStepMethodStep<Vector,
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

  ElastoDynamicsSolverStep(typename Base::EntitySet es,
                           std::shared_ptr<std::vector<int>> physical,
                           const Dune::ParameterTree& rootparams)
    : MaterialDependantStepBase<Vector>(es, physical, rootparams.sub("material"))
    , WrapperStep<Vector, OneStepMethodStep<Vector, SpatialLocalOperator, TemporalLocalOperator>>(std::make_shared<OneStepMethodStep<Vector, SpatialLocalOperator, TemporalLocalOperator>>())
    , params(rootparams)
  {}

  template<typename... FUNCS>
  ElastoDynamicsSolverStep(typename Base::EntitySet es,
                           std::shared_ptr<std::vector<int>> physical,
                           const Dune::ParameterTree& rootparams,
                           TimeCapsule<double>& tc,
                           const FUNCS&... funcs)
      : MaterialDependantStepBase<Vector>(es, physical, rootparams.sub("material"))
      , WrapperStep<Vector, OneStepMethodStep<Vector, SpatialLocalOperator, TemporalLocalOperator>>(std::make_shared<VariableBoundaryOneStepMethodStep<Vector, SpatialLocalOperator, TemporalLocalOperator>>(tc, funcs...))
      , params(rootparams)
  {}

  virtual ~ElastoDynamicsSolverStep() {}

  //TODO: Get rid of this which is necessary because of multiple inheritance
  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    this->step->update_parameter(name, param);
  }

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    auto& gfs = vector->gridFunctionSpace();
    this->step->set_spatial_localoperator(std::make_shared<SpatialLocalOperator>(gfs, gfs, params, this->material));
    this->step->set_temporal_localoperator(std::make_shared<TemporalLocalOperator>(gfs, gfs, params));
    this->step->pre(vector, cc);
  }

  //TODO: Get rid of this which is necessary because of multiple inheritance
  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    this->step->apply(vector, cc);
  }

  private:
  Dune::ParameterTree params;
};

#endif
