#ifndef DUNE_STRUCTURES_SOLVERSTEPS_REINFORCED_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_REINFORCED_HH

#include<dune/structures/eulerbernoulli.hh>
#include<dune/structures/solversteps/elasticity.hh>
#include<dune/structures/solversteps/traits.hh>


template<std::size_t num_fibres, typename... V>
class FibreReinforcedElasticitySolverStep
  : public WrapperStep<NewtonSolverTransitionStep<FibreReinforcedBulkOperator<typename SimpleStepTraits<V...>::GridFunctionSpace, SimpleStepTraits<V...>::dim>, V...>, V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  static constexpr int dim = Traits::dim;
  using LocalOperator = FibreReinforcedBulkOperator<typename Traits::GridFunctionSpace, dim>;

  FibreReinforcedElasticitySolverStep(const Dune::ParameterTree& rootparams, const Dune::ParameterTree& params)
    : WrapperStep<NewtonSolverTransitionStep<LocalOperator, V...>, V...>(std::make_shared<NewtonSolverTransitionStep<LocalOperator, V...>>(params))
    , rootparams(rootparams)
    , params(params)
  {}

  virtual ~FibreReinforcedElasticitySolverStep() {}

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
  {
    if (name == "material")
    {
      material = std::get<std::shared_ptr<typename Traits::Material>>(param);

      // I strongly dislike this lifetime dependencies, but Update parameter is available
      // very early right now - the LocalOperator might not even be constructed.
      if (lop != nullptr)
        lop->setMaterial(material);
    }

    this->step->update_parameter(name, param);
  }

  virtual void pre(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc) override
  {
    // Create the local operator for the bulk problem
    auto gfs = vector->gridFunctionSpaceStorage();
    lop = std::make_shared<LocalOperator>(gfs, rootparams, params, material);
    force = std::make_shared<typename Traits::Vector>(*gfs, 0.0);
    traction = std::make_shared<typename Traits::Vector>(*gfs, 0.0);
    lop->setCoefficientForce(gfs, force);
    lop->setCoefficientTraction(gfs, traction);

    this->step->set_localoperator(lop);
    this->step->pre(vector, cc);
  }

  virtual void apply(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc) override
  {
    auto gfs = vector->gridFunctionSpaceStorage();

    // Interpolate the force vector
    auto force_gf = makeGridFunctionTreeFromCallables(*gfs, get_callable_array<double(typename Traits::Entity, typename Traits::GlobalCoordinate), V...>(*(this->solver), params.get<std::string>("force", "0.0")));
    Dune::PDELab::interpolate(force_gf, *gfs, *force);

    // Interpolate the traction vector
    auto traction_gf = makeGridFunctionTreeFromCallables(*gfs, get_callable_array<double(typename Traits::Entity, typename Traits::GlobalCoordinate), V...>(*(this->solver), params.get<std::string>("traction", "0.0")));
    Dune::PDELab::interpolate(traction_gf, *gfs, *traction);

    this->step->apply(vector, cc);
  }

  private:
  std::shared_ptr<LocalOperator> lop;
  std::shared_ptr<typename Traits::Material> material;
  std::shared_ptr<typename Traits::Vector> force;
  std::shared_ptr<typename Traits::Vector> traction;
  Dune::ParameterTree rootparams;
  Dune::ParameterTree params;
};


#endif
