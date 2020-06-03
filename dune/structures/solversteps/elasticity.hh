#ifndef DUNE_STRUCTURES_SOLVERSTEPS_ELASTICITY_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_ELASTICITY_HH

#include<dune/structures/muparser.hh>
#include<dune/structures/elasticity.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/instationary.hh>
#include<dune/structures/solversteps/material.hh>
#include<dune/structures/solversteps/newton.hh>
#include<dune/structures/solversteps/traits.hh>


template<std::size_t i, typename... V>
class ElasticityOperatorStep
  : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = VectorStepTraits<i, V...>;

  static constexpr int dim = Traits::dim;
  using BaseOperator = AbstractLocalOperatorInterface<typename Traits::GridFunctionSpace>;
  using FGFS = typename std::tuple_element<i + 1, std::tuple<V...>>::type::GridFunctionSpace;
  using TGFS = typename std::tuple_element<i + 2, std::tuple<V...>>::type::GridFunctionSpace;
  using LocalOperator = typename OperatorSwitch<typename Traits::GridFunctionSpace, dim, FGFS, TGFS>::Elasticity;

  ElasticityOperatorStep(const Dune::ParameterTree& params)
    : params(params)
  {}

  virtual ~ElasticityOperatorStep() {}

  virtual void set_solver(std::shared_ptr<typename Traits::Solver> solver_) override
  {
    this->solver = solver_;

    // Construct the local operator...
    auto vector = this->solver->template getVector<i>();
    auto gfs = vector->gridFunctionSpaceStorage();
    auto material = this->solver->template param<std::shared_ptr<typename Traits::Material>>("material");
    auto lop = std::make_shared<LocalOperator>(*gfs, *gfs, params, material);

    auto force = this->solver->template getVector<i + 1>();
    lop->setCoefficientForce(force->gridFunctionSpaceStorage(), force);

    auto traction = this->solver->template getVector<i + 2>();
    lop->setCoefficientTraction(traction->gridFunctionSpaceStorage(), traction);

    // ... and register it in the parameter system
    this->solver->template introduce_parameter<std::shared_ptr<BaseOperator>>("elasticity_operator", lop);
  }

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
  {
    if (name == "material")
    {
      auto material = std::get<std::shared_ptr<typename Traits::Material>>(param);
      this->solver->template param<std::shared_ptr<BaseOperator>>("elasticity_operator")->setMaterial(material);
    }
  }

  private:
  Dune::ParameterTree params;
};


template<std::size_t i, typename... V>
class ElasticityMassOperatorStep
  : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = VectorStepTraits<i, V...>;

  using BaseOperator = AbstractLocalOperatorInterface<typename Traits::GridFunctionSpace>;
  using FGFS = typename std::tuple_element<i + 1, std::tuple<V...>>::type::GridFunctionSpace;
  using TGFS = typename std::tuple_element<i + 2, std::tuple<V...>>::type::GridFunctionSpace;
  using TemporalLocalOperator = typename OperatorSwitch<typename Traits::GridFunctionSpace,
                                                        Traits::dim, FGFS, TGFS>::Mass;

  ElasticityMassOperatorStep(const Dune::ParameterTree& params)
    : params(params)
  {}

  virtual ~ElasticityMassOperatorStep() override {}

  virtual void set_solver(std::shared_ptr<typename Traits::Solver> solver_) override
  {
    this->solver = solver_;

    // Construct the local operator...
    auto vector = this->solver->template getVector<i>();
    auto gfs = vector->gridFunctionSpaceStorage();
    auto lop = std::make_shared<TemporalLocalOperator>(*gfs, *gfs, params);

    // ... and register it in the parameter system
    this->solver->template introduce_parameter<std::shared_ptr<BaseOperator>>("elasticity_mass_operator", lop);
  }

  private:
  std::shared_ptr<typename Traits::Material> material;
  Dune::ParameterTree params;
};


#endif
