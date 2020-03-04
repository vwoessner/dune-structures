#ifndef DUNE_STRUCTURES_SOLVERSTEPS_REINFORCED_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_REINFORCED_HH

#include<dune/structures/eulerbernoulli.hh>
#include<dune/structures/solversteps/elasticity.hh>
#include<dune/structures/solversteps/traits.hh>


template<std::size_t i, typename... V>
class FibreReinforcedElasticitySolverStep
  : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = VectorStepTraits<i, V...>;

  static constexpr int dim = Traits::dim;
  using BaseOperator = AbstractLocalOperatorInterface<typename Traits::GridFunctionSpace>;
  using LocalOperator = FibreReinforcedBulkOperator<typename Traits::GridFunctionSpace, dim>;

  FibreReinforcedElasticitySolverStep(const Dune::ParameterTree& rootparams, const Dune::ParameterTree& params)
    : rootparams(rootparams)
    , params(params)
  {}

  virtual ~FibreReinforcedElasticitySolverStep() {}

  virtual void set_solver(std::shared_ptr<typename Traits::Solver> solver_)
  {
    this->solver = solver_;

    // Construct the local operator...
    auto vector = this->solver->template getVector<i>();
    auto gfs = vector->gridFunctionSpaceStorage();
    auto material = this->solver->template param<std::shared_ptr<typename Traits::Material>>("material");
    auto lop = std::make_shared<LocalOperator>(gfs, rootparams, params, material);
    auto force = std::make_shared<typename Traits::Vector>(*gfs, 0.0);
    auto traction = std::make_shared<typename Traits::Vector>(*gfs, 0.0);

    // Interpolate the force vector
    auto force_gf = makeGridFunctionTreeFromCallables(*gfs, get_callable_array<double(typename Traits::Entity, typename Traits::GlobalCoordinate), V...>(*(this->solver), params.get<std::string>("force", "0.0")));
    Dune::PDELab::interpolate(force_gf, *gfs, *force);

    // Interpolate the traction vector
    auto traction_gf = makeGridFunctionTreeFromCallables(*gfs, get_callable_array<double(typename Traits::Entity, typename Traits::GlobalCoordinate), V...>(*(this->solver), params.get<std::string>("traction", "0.0")));
    Dune::PDELab::interpolate(traction_gf, *gfs, *traction);

    lop->setCoefficientForce(gfs, force);
    lop->setCoefficientTraction(gfs, traction);

    // ... and register it in the parameter system
    this->solver->template introduce_parameter<std::shared_ptr<BaseOperator>>("fibre_operator", lop);
  }

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
  {
    if (name == "material")
    {
      auto material = std::get<std::shared_ptr<typename Traits::Material>>(param);
      this->solver->template param<std::shared_ptr<BaseOperator>>("fibre_operator")->setMaterial(material);
    }
  }

  private:
  Dune::ParameterTree rootparams;
  Dune::ParameterTree params;
};

#endif
