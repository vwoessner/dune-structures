#ifndef DUNE_STRUCTURES_SOLVERSTEPS_ELASTICITY_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_ELASTICITY_HH

#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/material.hh>
#include<dune/structures/solversteps/newton.hh>

#include"operators/elasticity_operator.hh"
#include"operators/elastodynamics_spatial_operator.hh"
#include"operators/elastodynamics_temporal_operator.hh"


template<typename Vector>
class ElasticitySolverStep
  : public MaterialDependantStepBase<Vector>
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
    , newton_step(std::make_shared<NewtonSolverTransitionStep<Vector, LocalOperator>>())
    , params(rootparams)
  {}

  virtual ~ElasticitySolverStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    // Parse the material
    auto& gfs = vector->gridFunctionSpace();

    // Instantiate the local operator
    newton_step->set_localoperator(std::make_shared<LocalOperator>(gfs, gfs, params, this->material));

    // Call the base class
    newton_step->apply(vector, cc);
  }

  private:
  std::shared_ptr<NewtonSolverTransitionStep<Vector, LocalOperator>> newton_step;
  Dune::ParameterTree params;
};

#endif
