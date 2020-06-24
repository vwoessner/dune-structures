#ifndef DUNE_STRUCTURES_ELASTICITY_HH
#define DUNE_STRUCTURES_ELASTICITY_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/blocks/enableif.hh>
#include<dune/pdelab.hh>

#include<iostream>
#include<memory>
#include<tuple>

// The 3D operators
#include"operators/elasticity_operator.hh"
#include"operators/quasistatic_mass_operator.hh"
#include"operators/elasticity_p2_operator.hh"
#include"operators/quasistatic_mass_p2_operator.hh"

// The 2D operators
#include"operators/elasticity_2d_operator.hh"
#include"operators/quasistatic_mass_2d_operator.hh"
#include"operators/elasticity_2d_p2_operator.hh"
#include"operators/quasistatic_mass_2d_p2_operator.hh"


template<typename GFS, int dim, typename FEM, typename FGFS, typename TGFS>
struct OperatorSwitchImpl
{};

template<typename GFS, int dim, typename FGFS, typename TGFS>
using OperatorSwitch = OperatorSwitchImpl<GFS, dim, typename GFS::template Child<0>::Type::Traits::FiniteElementMap, FGFS, TGFS>;

template<typename GFS, typename ES, typename R, typename FGFS, typename TGFS>
struct OperatorSwitchImpl<GFS, 3, Dune::PDELab::PkLocalFiniteElementMap<ES, double, R, 1>, FGFS, TGFS>
{
  using Elasticity = ElasticityOperator<GFS, GFS, FGFS, TGFS>;
  using Mass = QuasiStaticMassOperator<GFS, GFS>;
};

template<typename GFS, typename ES, typename R, typename FGFS, typename TGFS>
struct OperatorSwitchImpl<GFS, 3, Dune::PDELab::PkLocalFiniteElementMap<ES, double, R, 2>, FGFS, TGFS>
{
  using Elasticity = ElasticityP2Operator<GFS, GFS, FGFS, TGFS>;
  using Mass = QuasiStaticMassP2Operator<GFS, GFS>;
};

template<typename GFS, typename ES, typename R, typename FGFS, typename TGFS>
struct OperatorSwitchImpl<GFS, 2, Dune::PDELab::PkLocalFiniteElementMap<ES, double, R, 1>, FGFS, TGFS>
{
  using Elasticity = Elasticity2DOperator<GFS, GFS, FGFS, TGFS>;
  using Mass = QuasiStaticMass2DOperator<GFS, GFS>;
};

template<typename GFS, typename ES, typename R, typename FGFS, typename TGFS>
struct OperatorSwitchImpl<GFS, 2, Dune::PDELab::PkLocalFiniteElementMap<ES, double, R, 2>, FGFS, TGFS>
{
  using Elasticity = Elasticity2DP2Operator<GFS, GFS, FGFS, TGFS>;
  using Mass = QuasiStaticMass2DP2Operator<GFS, GFS>;
};


template<typename P, typename V, std::size_t i, typename enable = Dune::BlockLab::disabled>
class ElasticityOperatorBlock
  : public Dune::BlockLab::DisabledBlock<P, V, i>
{
  public:
  template<typename Context>
  ElasticityOperatorBlock(Context& ctx, const Dune::ParameterTree& config)
    : Dune::BlockLab::DisabledBlock<P, V, i>(ctx, config)
  {}
};


template<typename P, typename V, std::size_t i>
class ElasticityOperatorBlock<P, V, i, Dune::BlockLab::enableBlock<Dune::BlockLab::accessesAdditionalVectors<P, V, i, 2>()>>
  : public Dune::BlockLab::BlockBase<P, V, i>
{
  public:
  using Traits = Dune::BlockLab::BlockTraits<P, V, i>;
  using Material = std::shared_ptr<ElasticMaterialBase<typename Traits::EntitySet, double>>;

  static constexpr int dim = Traits::dim;
  using BaseOperator = Dune::BlockLab::AbstractLocalOperatorInterface<typename Traits::GridFunctionSpace>;
  using FGFS = typename std::tuple_element<i + 1, V>::type::GridFunctionSpace;
  using TGFS = typename std::tuple_element<i + 2, V>::type::GridFunctionSpace;
  using LocalOperator = typename OperatorSwitch<typename Traits::GridFunctionSpace, dim, FGFS, TGFS>::Elasticity;

  template<typename Context>
  ElasticityOperatorBlock(Context&, const Dune::ParameterTree& config)
    : params(config)
  {}

  virtual ~ElasticityOperatorBlock() = default;

  virtual void setup() override
  {
    // Construct the local operator...
    auto vector = this->solver->template getVector<i>();
    auto gfs = vector->gridFunctionSpaceStorage();
    auto material = this->solver->template param<Material>("material");
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
      auto material = std::get<Material>(param);
      auto lop_pointer = this->solver->template param<std::shared_ptr<BaseOperator>>("elasticity_operator").get();
      dynamic_cast<LocalOperator*>(lop_pointer)->setMaterial(material);
    }
  }

  private:
  Dune::ParameterTree params;
};



template<typename P, typename V, std::size_t i, typename enable = Dune::BlockLab::disabled>
class ElasticityMassOperatorBlock
  : public Dune::BlockLab::DisabledBlock<P, V, i>
{
  public:
  template<typename Context>
  ElasticityMassOperatorBlock(Context& ctx, const Dune::ParameterTree& config)
    : Dune::BlockLab::DisabledBlock<P, V, i>(ctx, config)
  {}
};

template<typename P, typename V, std::size_t i>
class ElasticityMassOperatorBlock<P, V, i, Dune::BlockLab::enableBlock<Dune::BlockLab::accessesAdditionalVectors<P, V, i, 2>()>>
  : public Dune::BlockLab::BlockBase<P, V, i>
{
  public:
  using Traits = Dune::BlockLab::BlockTraits<P, V, i>;

  using BaseOperator = Dune::BlockLab::AbstractLocalOperatorInterface<typename Traits::GridFunctionSpace>;
  using FGFS = typename std::tuple_element<i + 1, V>::type::GridFunctionSpace;
  using TGFS = typename std::tuple_element<i + 2, V>::type::GridFunctionSpace;
  using TemporalLocalOperator = typename OperatorSwitch<typename Traits::GridFunctionSpace,
                                                        Traits::dim, FGFS, TGFS>::Mass;

  template<typename Context>
  ElasticityMassOperatorBlock(Context&, const Dune::ParameterTree& config)
    : params(config)
  {}

  virtual ~ElasticityMassOperatorBlock() = default;

  virtual void setup() override
  {
    // Construct the local operator...
    auto vector = this->solver->template getVector<i>();
    auto gfs = vector->gridFunctionSpaceStorage();
    auto lop = std::make_shared<TemporalLocalOperator>(*gfs, *gfs, params);

    // ... and register it in the parameter system
    this->solver->template introduce_parameter<std::shared_ptr<BaseOperator>>("elasticity_mass_operator", lop);
  }

  private:
  Dune::ParameterTree params;
};

#endif
