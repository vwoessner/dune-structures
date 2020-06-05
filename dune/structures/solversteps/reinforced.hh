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
  using FGFS = typename std::tuple_element<i + 1, std::tuple<V...>>::type::GridFunctionSpace;
  using TGFS = typename std::tuple_element<i + 2, std::tuple<V...>>::type::GridFunctionSpace;
  using BaseOperator = AbstractLocalOperatorInterface<typename Traits::GridFunctionSpace>;
  using LocalOperator = FibreReinforcedBulkOperator<typename Traits::GridFunctionSpace, FGFS, TGFS, dim>;

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

    auto force = this->solver->template getVector<i + 1>();
    lop->setCoefficientForce(force->gridFunctionSpaceStorage(), force);

    auto traction = this->solver->template getVector<i + 2>();
    lop->setCoefficientTraction(traction->gridFunctionSpaceStorage(), traction);

    // ... and register it in the parameter system
    this->solver->template introduce_parameter<std::shared_ptr<BaseOperator>>("fibre_operator", lop);
  }

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
  {
    if (name == "material")
    {
      auto material = std::get<std::shared_ptr<typename Traits::Material>>(param);
      auto lop_pointer = this->solver->template param<std::shared_ptr<BaseOperator>>("fibre_operator").get();
      dynamic_cast<LocalOperator*>(lop_pointer)->setMaterial(material);
    }
    if (name == "adapted")
    {
      auto lop_pointer = this->solver->template param<std::shared_ptr<BaseOperator>>("fibre_operator").get();
      dynamic_cast<LocalOperator*>(lop_pointer)->compute_grid_intersection();
    }
  }

  private:
  Dune::ParameterTree rootparams;
  Dune::ParameterTree params;
};

#endif
