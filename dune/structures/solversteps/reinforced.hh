#ifndef DUNE_STRUCTURES_SOLVERSTEPS_REINFORCED_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_REINFORCED_HH

#include<dune/structures/eulerbernoulli.hh>
#include<dune/structures/solversteps/elasticity.hh>


template<typename Vector, std::size_t num_fibres>
class FibreReinforcedElasticitySolverStep
  : public WrapperStep<Vector,
                       NewtonSolverTransitionStep<Vector, FibreReinforcedBulkOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace, TransitionSolverStepBase<Vector>::dim>>>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;
  static constexpr int dim = Base::dim;
  using LocalOperator = FibreReinforcedBulkOperator<typename Base::GridFunctionSpace, dim>;

  FibreReinforcedElasticitySolverStep(const Dune::ParameterTree& rootparams, const Dune::ParameterTree& params)
    : WrapperStep<Vector, NewtonSolverTransitionStep<Vector, LocalOperator>>(std::make_shared<NewtonSolverTransitionStep<Vector, LocalOperator>>())
    , rootparams(rootparams)
    , params(params)
  {}

  virtual ~FibreReinforcedElasticitySolverStep() {}

  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    if (name == "material")
    {
      material = std::get<std::shared_ptr<typename Base::Material>>(param);

      // I strongly dislike this lifetime dependencies, but Update parameter is available
      // very early right now - the LocalOperator might not even be constructed.
      if (lop != nullptr)
        lop->setMaterial(material);
    }

    this->step->update_parameter(name, param);
  }

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    // Create the local operator for the bulk problem
    auto gfs = vector->gridFunctionSpaceStorage();
    lop = std::make_shared<LocalOperator>(gfs, rootparams, params, material);
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
  std::shared_ptr<LocalOperator> lop;
  std::shared_ptr<typename Base::Material> material;
  std::shared_ptr<Vector> force;
  std::shared_ptr<Vector> traction;
  Dune::ParameterTree rootparams;
  Dune::ParameterTree params;
};


#endif
