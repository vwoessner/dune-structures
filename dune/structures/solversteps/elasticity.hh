#ifndef DUNE_STRUCTURES_SOLVERSTEPS_ELASTICITY_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_ELASTICITY_HH

#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/newton.hh>

#include"elasticity_operator.hh"


template<typename Vector>
class ElasticitySolverStep
  : public NewtonSolverTransitionStep<Vector,
                                      ElasticityOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                         typename TransitionSolverStepBase<Vector>::GridFunctionSpace>
                                      >
{
  public:
  using Base = TransitionSolverStepBase<Vector>;
  using LocalOperator = ElasticityOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                           typename TransitionSolverStepBase<Vector>::GridFunctionSpace>;

  ElasticitySolverStep(std::shared_ptr<std::vector<int>> physical,
                       const Dune::ParameterTree& params
                       )
    : physical(physical),
      params(params),
      material_params(params.sub("material"))
  {}

  virtual ~ElasticitySolverStep() {}

  virtual void apply(Vector& vector, typename Base::ConstraintsContainer& cc) override
  {
    // Parse the material
    auto gfs = vector.gridFunctionSpace();
    auto es = gfs.entitySet();
    auto material = parse_material<double>(es, physical, material_params);

    // Instantiate the local operator
    this->localoperator = std::make_shared<LocalOperator>(gfs, gfs, params, material);

    // Call the base class
    this->NewtonSolverTransitionStep<Vector, LocalOperator>::apply(vector, cc);
  }

  private:
  std::shared_ptr<std::vector<int>> physical;
  const Dune::ParameterTree& params;
  Dune::ParameterTree material_params;
};

#endif
