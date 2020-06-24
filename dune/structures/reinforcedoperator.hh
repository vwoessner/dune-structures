#ifndef DUNE_STRUCTURES_SOLVERSTEPS_REINFORCED_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_REINFORCED_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/blocks/enableif.hh>
#include<dune/blocklab/operators/virtualinterface.hh>
#include<dune/common/parametertree.hh>
#include<dune/structures/eulerbernoulli.hh>


template<typename P, typename V, std::size_t i, typename enable = Dune::BlockLab::disabled>
class FibreReinforcedElasticityOperatorBlock
  : public Dune::BlockLab::DisabledBlock<P, V, i>
{
  public:
  template<typename Context>
  FibreReinforcedElasticityOperatorBlock(Context& ctx, const Dune::ParameterTree& config)
    : Dune::BlockLab::DisabledBlock<P, V, i>(ctx, config)
  {}
};


template<typename P, typename V, std::size_t i>
class FibreReinforcedElasticityOperatorBlock<P, V, i, Dune::BlockLab::enableBlock<Dune::BlockLab::accessesAdditionalVectors<P, V, i, 2>() && Dune::BlockLab::is2D<P, V, i>()>>
  : public Dune::BlockLab::BlockBase<P, V, i>
{
  public:
  using Traits = Dune::BlockLab::BlockTraits<P, V, i>;
  using Material = std::shared_ptr<ElasticMaterialBase<typename Traits::EntitySet, double>>;

  static constexpr int dim = Traits::dim;
  using FGFS = typename std::tuple_element<i + 1, V>::type::GridFunctionSpace;
  using TGFS = typename std::tuple_element<i + 2, V>::type::GridFunctionSpace;
  using BaseOperator = Dune::BlockLab::AbstractLocalOperatorInterface<typename Traits::GridFunctionSpace>;
  using LocalOperator = FibreReinforcedBulkOperator<typename Traits::GridFunctionSpace, FGFS, TGFS, dim>;

  template<typename Context>
  FibreReinforcedElasticityOperatorBlock(Context& ctx, const Dune::ParameterTree& config)
    : rootparams(ctx.getRootConfig())
    , params(config)
  {}

  virtual ~FibreReinforcedElasticityOperatorBlock() = default;

  virtual void setup() override
  {
    // Construct the local operator...
    auto vector = this->solver->template getVector<i>();
    auto gfs = vector->gridFunctionSpaceStorage();
    auto material = this->solver->template param<Material>("material");
    auto lop = std::make_shared<LocalOperator>(gfs, rootparams, params, material);

    auto force = this->solver->template getVector<i + 1>();
    lop->setCoefficientForce(force->gridFunctionSpaceStorage(), force);

    auto traction = this->solver->template getVector<i + 2>();
    lop->setCoefficientTraction(traction->gridFunctionSpaceStorage(), traction);

    // ... and register it in the parameter system
    this->solver->template introduce_parameter<std::shared_ptr<BaseOperator>>("reinforced_operator", lop);
  }

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
  {
    if (name == "material")
    {
      auto material = std::get<Material>(param);
      auto lop_pointer = this->solver->template param<std::shared_ptr<BaseOperator>>("reinforced_operator").get();
      dynamic_cast<LocalOperator*>(lop_pointer)->setMaterial(material);
    }
    if (name == "adapted")
    {
      auto lop_pointer = this->solver->template param<std::shared_ptr<BaseOperator>>("reinforced_operator").get();
      dynamic_cast<LocalOperator*>(lop_pointer)->compute_grid_intersection();
    }
  }

  private:
  Dune::ParameterTree rootparams;
  Dune::ParameterTree params;
};

#endif
