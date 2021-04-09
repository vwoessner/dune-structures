#ifndef DUNE_STRUCTURES_SOLVERSTEPS_REINFORCED_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_REINFORCED_HH

#include <dune/blocklab/blocks/blockbase.hh>
#include <dune/blocklab/blocks/enableif.hh>
#include <dune/blocklab/operators/virtualinterface.hh>
#include <dune/common/parametertree.hh>
#include <dune/structures/eulerbernoulli.hh>

template<typename P,
         typename V,
         std::size_t i,
         typename enable = Dune::BlockLab::disabled>
class FibreReinforcedElasticityOperatorBlock
  : public Dune::BlockLab::DisabledBlock<P, V, i>
{
public:
  template<typename Context>
  FibreReinforcedElasticityOperatorBlock(Context& ctx, const YAML::Node& config)
    : Dune::BlockLab::DisabledBlock<P, V, i>(ctx, config)
  {
  }
};

template<typename P, typename V, std::size_t i>
class FibreReinforcedElasticityOperatorBlock<
  P,
  V,
  i,
  Dune::BlockLab::enableBlock<
    Dune::BlockLab::accessesAdditionalVectors<P, V, i, 2>()
    && Dune::BlockLab::is2D<P, V, i>()>>
  : public Dune::BlockLab::BlockBase<P, V, i>
{
public:
  using Traits = Dune::BlockLab::BlockTraits<P, V, i>;
  using Material =
    std::shared_ptr<ElasticMaterialBase<typename Traits::EntitySet, double>>;

  static constexpr int dim = Traits::dim;
  using FGFS = typename std::tuple_element<i + 1, V>::type::GridFunctionSpace;
  using TGFS = typename std::tuple_element<i + 2, V>::type::GridFunctionSpace;
  using BaseOperator = Dune::BlockLab::AbstractLocalOperatorInterface<
    typename Traits::GridFunctionSpace>;
  using LocalOperator =
    FibreReinforcedBulkOperator<typename Traits::GridFunctionSpace,
                                FGFS,
                                TGFS,
                                dim>;

  template<typename Context>
  FibreReinforcedElasticityOperatorBlock(Context& ctx, const YAML::Node& config)
    : Dune::BlockLab::BlockBase<P, V, i>(ctx, config)
    , params(config)
  {
  }

  virtual ~FibreReinforcedElasticityOperatorBlock() = default;

  virtual void setup() override
  {
    // Construct the local operator...
    auto vector = this->solver->template getVector<i>();
    auto gfs = vector->gridFunctionSpaceStorage();
    auto material = this->solver->template param<Material>("material");
    auto lop = std::make_shared<LocalOperator>(gfs, params, material);

    auto force = this->solver->template getVector<i + 1>();
    lop->setCoefficientForce(force->gridFunctionSpaceStorage(), force);

    auto traction = this->solver->template getVector<i + 2>();
    lop->setCoefficientTraction(traction->gridFunctionSpaceStorage(), traction);

    // ... and register it in the parameter system
    this->solver->template introduce_parameter<std::shared_ptr<BaseOperator>>(
      this->getBlockName(), lop);
  }

  virtual void update_parameter(std::string name,
                                typename Traits::Parameter param) override
  {
    if (name == "material")
    {
      auto material = std::get<Material>(param);
      auto lop_pointer =
        this->solver
          ->template param<std::shared_ptr<BaseOperator>>(this->getBlockName())
          .get();
      dynamic_cast<LocalOperator*>(lop_pointer)->setMaterial(material);
    }
    if (name == "adapted")
    {
      auto lop_pointer =
        this->solver
          ->template param<std::shared_ptr<BaseOperator>>(this->getBlockName())
          .get();
      dynamic_cast<LocalOperator*>(lop_pointer)->compute_grid_intersection();
    }
  }

  static std::vector<std::string> blockData()
  {
    auto data = Dune::BlockLab::BlockBase<P, V, i>::blockData();
    data.push_back("title: Fibre-reinforced Elasticity Operator         \n"
                   "category: operators                                 \n"
                   "schema:                                             \n"
                   "  stabilization_parameter:                          \n"
                   "    type: float                                     \n"
                   "    default: 1.0                                    \n"
                   "    meta:                                           \n"
                   "      title: DG Stabilization Parameter             \n"
                   "  fibres:                                           \n"
                   "    type: list                                      \n"
                   "    schema:                                         \n"
                   "      type: dict                                    \n"
                   "      schema:                                       \n"
                   "        start:                                      \n"
                   "          type: list                                \n"
                   "          maxlength: 2                              \n"
                   "          minlength: 2                              \n"
                   "          schema:                                   \n"
                   "            type: float                             \n"
                   "            default: 0                              \n"
                   "          meta:                                     \n"
                   "            title: Start point                      \n"
                   "        end:                                        \n"
                   "          type: list                                \n"
                   "          maxlength: 2                              \n"
                   "          minlength: 2                              \n"
                   "          schema:                                   \n"
                   "            type: float                             \n"
                   "            default: 1                              \n"
                   "          meta:                                     \n"
                   "            title: End point                        \n"
                   "        radius:                                     \n"
                   "          type: float                               \n"
                   "          default: 0.1                              \n"
                   "          meta:                                     \n"
                   "            title: Radius                           \n"
                   "        youngs_modulus:                             \n"
                   "          type: float                               \n"
                   "          default: 1e5                              \n"
                   "          meta:                                     \n"
                   "            title: Young's modulus                  \n"
                   "    meta:                                           \n"
                   "      title: Fibres                                 \n");
    return data;
  }

private:
  YAML::Node params;
};

#endif
