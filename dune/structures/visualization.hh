#ifndef DUNE_STRUCTURES_VISUALIZATION_HH
#define DUNE_STRUCTURES_VISUALIZATION_HH

#include <dune/blocklab/blocks/blockbase.hh>
#include <dune/blocklab/blocks/enableif.hh>
#include <dune/blocklab/blocks/visualization.hh>
#include <dune/blocklab/utilities/yaml.hh>
#include <dune/structures/material.hh>
#include <dune/structures/stress.hh>
#include <dune/structures/vonmises.hh>

#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/test/l2norm.hh>

#include <filesystem>
#include <memory>
#include <variant>
#include <vector>

template<typename P, typename V>
class PhysicalEntityVisualizationBlock : public Dune::BlockLab::BlockBase<P, V>
{
public:
  using Traits = Dune::BlockLab::BlockTraits<P, V>;

  template<typename Context>
  PhysicalEntityVisualizationBlock(Context& ctx, const YAML::Node& config)
    : Dune::BlockLab::BlockBase<P, V>(ctx, config)
  {
  }

  virtual ~PhysicalEntityVisualizationBlock() = default;

  struct PhysicalEntityFunction
    : public Dune::VTKFunction<typename Traits::GridView>
  {
    using Base = Dune::VTKFunction<typename Traits::GridView>;
    using Entity = typename Base::Entity;

    PhysicalEntityFunction(typename Traits::GridView gv,
                           std::shared_ptr<std::vector<int>> physical)
      : is(gv.indexSet())
      , physical(physical)
    {
    }

    virtual ~PhysicalEntityFunction() = default;

    virtual int ncomps() const override { return 1; }

    virtual double evaluate(
      int comp,
      const typename Traits::Entity& e,
      const typename Traits::LocalCoordinate& xi) const override
    {
      return static_cast<double>((*physical)[is.index(e)]);
    }

    virtual std::string name() const override
    {
      return "GMSH Physical Entity Information";
    }

    const typename Traits::GridView::IndexSet& is;
    std::shared_ptr<std::vector<int>> physical;
  };

  virtual void setup() override
  {
    auto gv =
      this->solver->template getVector<0>()->gridFunctionSpace().gridView();
    auto physical =
      this->solver->template param<std::shared_ptr<std::vector<int>>>(
        "physical");
    auto function = std::make_shared<PhysicalEntityFunction>(gv, physical);
    std::dynamic_pointer_cast<Dune::BlockLab::VisualizationBlock<P, V>>(
      this->parent)
      ->add_celldata(function);
  }

  static std::vector<std::string> blockData()
  {
    auto data = Dune::BlockLab::BlockBase<P, V>::blockData();
    data.push_back("title: Physical entity information visualization    \n"
                   "category: visualization                             \n");
    return data;
  }
};

template<typename P,
         typename V,
         std::size_t i,
         typename enabled = Dune::BlockLab::disabled>
class VonMisesStressVisualizationBlock
  : public Dune::BlockLab::DisabledBlock<P, V, i>
{
public:
  template<typename Context>
  VonMisesStressVisualizationBlock(Context& ctx, const YAML::Node& config)
    : Dune::BlockLab::DisabledBlock<P, V, i>(ctx, config)
  {
  }
};

template<typename P, typename V, std::size_t i>
class VonMisesStressVisualizationBlock<
  P,
  V,
  i,
  Dune::BlockLab::enableBlock<Dune::BlockLab::isDimPower<P, V, i>()>>
  : public Dune::BlockLab::BlockBase<P, V, i>
{
public:
  using Traits = Dune::BlockLab::BlockTraits<P, V, i>;
  using Material =
    std::shared_ptr<ElasticMaterialBase<typename Traits::EntitySet, double>>;

  template<typename Context>
  VonMisesStressVisualizationBlock(Context& ctx, const YAML::Node& config)
    : Dune::BlockLab::BlockBase<P, V, i>(ctx, config)
    , continuous(config["continuous"].as<bool>())
    , report_l2(config["report_l2"].as<bool>())
    , report_inf(config["report_inf"].as<bool>())
  {
  }

  virtual ~VonMisesStressVisualizationBlock() = default;

  virtual void update_parameter(std::string name,
                                typename Traits::Parameter param) override
  {
    if (name == "material")
      material = std::get<Material>(param);
  }

  virtual void apply() override
  {
    auto vector = this->solver->template getVector<i>();
    auto es = vector->gridFunctionSpace().entitySet();
    auto gfs = vector->gridFunctionSpaceStorage();
    gfs->ordering();

    // A grid function for the stress
    VonMisesStressGridFunction<typename Traits::Vector, Traits::dim> stress(
      *vector, material);

    // Interpolate the stress into a grid function
    if (continuous)
    {
      using P1FEM =
        Dune::PDELab::PkLocalFiniteElementMap<typename Traits::EntitySet,
                                              double,
                                              typename Traits::Range,
                                              1>;
      auto p1fem = std::make_shared<P1FEM>(es);
      using P1GFS =
        Dune::PDELab::GridFunctionSpace<typename Traits::EntitySet,
                                        P1FEM,
                                        Dune::PDELab::NoConstraints,
                                        typename Traits::VectorBackend>;
      auto p1gfs = std::make_shared<P1GFS>(es, p1fem);
      using StressVector =
        Dune::PDELab::Backend::Vector<P1GFS, typename Traits::ctype>;
      auto stress_container = std::make_shared<StressVector>(p1gfs);

      Dune::PDELab::interpolate(stress, *p1gfs, *stress_container);
      std::dynamic_pointer_cast<Dune::BlockLab::VisualizationBlock<P, V>>(
        this->parent)
        ->add_dataset(stress_container, "vonmises");

      if (report_inf)
      {
        report_inf_to_cout(stress_container);
      }
    }
    else
    {
      using P0FEM =
        Dune::PDELab::P0LocalFiniteElementMap<typename Traits::ctype,
                                              typename Traits::Range,
                                              Traits::dim>;
      auto fem =
        std::make_shared<P0FEM>(Dune::GeometryTypes::simplex(Traits::dim));
      using P0GFS =
        Dune::PDELab::GridFunctionSpace<typename Traits::EntitySet,
                                        P0FEM,
                                        Dune::PDELab::NoConstraints,
                                        typename Traits::VectorBackend>;
      auto p0gfs = std::make_shared<P0GFS>(es, fem);
      p0gfs->setDataSetType(
        Dune::PDELab::GridFunctionOutputParameters::Output::cellData);
      using StressVector =
        Dune::PDELab::Backend::Vector<P0GFS, typename Traits::ctype>;
      auto stress_container = std::make_shared<StressVector>(p0gfs);

      Dune::PDELab::interpolate(stress, *p0gfs, *stress_container);
      std::dynamic_pointer_cast<Dune::BlockLab::VisualizationBlock<P, V>>(
        this->parent)
        ->add_dataset(stress_container, "vonmises");

      if (report_inf)
      {
        report_inf_to_cout(stress_container);
      }
    }

    if (report_l2)
    {
      std::cout << std::scientific << "stress l2: " << l2norm(stress, 2)
                << std::defaultfloat << std::endl;
    }
  }

  static std::vector<std::string> blockData()
  {
    auto data = Dune::BlockLab::BlockBase<P, V, i>::blockData();
    data.push_back("title: Von-Mises Stress Visualization               \n"
                   "category: visualization                             \n"
                   "schema:                                             \n"
                   "  continuous:                                       \n"
                   "    type: boolean                                   \n"
                   "    default: false                                  \n"
                   "    meta:                                           \n"
                   "      title: Use Continuous Interpolation           \n"
                   "  report_l2:                                        \n"
                   "    type: boolean                                   \n"
                   "    default: false                                  \n"
                   "    meta:                                           \n"
                   "      title: Report the global L2 norm of the stress\n"
                   "  report_inf:                                       \n"
                   "    type: boolean                                   \n"
                   "    default: false                                  \n"
                   "    meta:                                           \n"
                   "      title: Report the infinity norm of the stress \n");
    return data;
  }

private:
  Material material;
  bool continuous, report_l2, report_inf;

  /// Report the infinity norm of a backend vector to the terminal
  template<class BackendVector>
  void report_inf_to_cout(const std::shared_ptr<BackendVector> vec)
  {
    std::cout << std::scientific << "stress inf: " << vec->infinity_norm()
              << std::defaultfloat << std::endl;
  }
};

template<typename P, typename V>
class FibreDistanceVisualizationBlock : public Dune::BlockLab::BlockBase<P, V>
{
public:
  using Traits = Dune::BlockLab::BlockTraits<P, V>;

  template<typename Context>
  FibreDistanceVisualizationBlock(Context& ctx, const YAML::Node& config)
    : Dune::BlockLab::BlockBase<P, V>(ctx, config)
    , prestress(config, ctx.getRootConfig())
  {
  }

  virtual ~FibreDistanceVisualizationBlock() = default;

  virtual void setup() override
  {
    if constexpr (Traits::dim == 3)
    {
      auto vector = this->solver->template getVector<0>();
      auto es = vector->gridFunctionSpace().entitySet();
      using FEM =
        Dune::PDELab::PkLocalFiniteElementMap<typename Traits::EntitySet,
                                              double,
                                              typename Traits::Range,
                                              1>;
      auto fem = std::make_shared<FEM>(es);
      using GFS =
        Dune::PDELab::GridFunctionSpace<typename Traits::EntitySet,
                                        FEM,
                                        Dune::PDELab::NoConstraints,
                                        typename Traits::VectorBackend>;
      auto gfs = std::make_shared<GFS>(es, fem);
      gfs->name("fibredistance");
      using DistanceVector =
        Dune::PDELab::Backend::Vector<GFS, typename Traits::ctype>;
      auto container = std::make_shared<DistanceVector>(gfs);

      auto lambda = [this](const auto& e, const auto& x) {
        return this->prestress.distance_to_minimum(e, x);
      };
      auto gf =
        Dune::PDELab::makeGridFunctionFromCallable(es.gridView(), lambda);
      Dune::PDELab::interpolate(gf, *gfs, *container);
      std::dynamic_pointer_cast<Dune::BlockLab::VisualizationBlock<P, V>>(
        this->parent)
        ->add_dataset(container, "fibredistance");
    }
  }

  static std::vector<std::string> blockData()
  {
    auto data = Dune::BlockLab::BlockBase<P, V>::blockData();
    data.push_back("title: Distance to next fibre Visualization         \n"
                   "category: visualization                             \n");
    return data;
  }

private:
  CurvedFibrePrestress<typename Traits::GridView, double> prestress;
};

/// Visualize stress tensor eigenvectors
/** Disable this block if the vector it is applied to does not have the power
 *  of the dimension.
 */
template<typename P,
         typename V,
         std::size_t i,
         typename enabled = Dune::BlockLab::disabled>
class StressEVVisualizationBlock : public Dune::BlockLab::DisabledBlock<P, V, i>
{
public:
  template<typename Context>
  StressEVVisualizationBlock(Context& ctx, const YAML::Node& config)
    : Dune::BlockLab::DisabledBlock<P, V, i>(ctx, config)
  {
  }
};

/// Visualize stress tensor eigenvectors
template<typename P, typename V, std::size_t i>
class StressEVVisualizationBlock<
  P,
  V,
  i,
  Dune::BlockLab::enableBlock<Dune::BlockLab::isDimPower<P, V, i>()>>
  : public Dune::BlockLab::BlockBase<P, V, i>
{
public:
  using Traits = Dune::BlockLab::BlockTraits<P, V, i>;
  using Material =
    std::shared_ptr<ElasticMaterialBase<typename Traits::EntitySet, double>>;

  template<typename Context>
  StressEVVisualizationBlock(Context& ctx, const YAML::Node& config)
    : Dune::BlockLab::BlockBase<P, V, i>(ctx, config)
    , _index(config["index"].as<int>())
    , _continuous(config["continuous"].as<bool>())
  {
  }

  virtual ~StressEVVisualizationBlock() = default;

  virtual void update_parameter(std::string name,
                                typename Traits::Parameter param) override
  {
    if (name == "material")
      material = std::get<Material>(param);
  }

  virtual void apply() override
  {
    auto vector = this->solver->template getVector<i>();
    auto es = vector->gridFunctionSpace().entitySet();

    StressEVGridFunction<typename Traits::Vector, Traits::dim> stress(*vector,
                                                                      material);
    stress.set_index(_index);

    if (_continuous)
    {
      // P1 representation
      using FEM =
        Dune::PDELab::PkLocalFiniteElementMap<typename Traits::EntitySet,
                                              double,
                                              typename Traits::Range,
                                              1>;
      auto fem = std::make_shared<FEM>(es);

      using VGFS =
        Dune::PDELab::VectorGridFunctionSpace<typename Traits::EntitySet,
                                              FEM,
                                              Traits::dim,
                                              typename Traits::VectorBackend,
                                              typename Traits::VectorBackend,
                                              Dune::PDELab::NoConstraints>;
      auto vgfs = std::make_shared<VGFS>(es, fem);
      using StressVector =
        Dune::PDELab::Backend::Vector<VGFS, typename Traits::ctype>;
      auto stress_container = std::make_shared<StressVector>(vgfs);

      Dune::PDELab::interpolate(stress, *vgfs, *stress_container);
      std::dynamic_pointer_cast<Dune::BlockLab::VisualizationBlock<P, V>>(
        this->parent)
        ->add_dataset(stress_container, "stress_ev_" + std::to_string(_index));
    }
    else
    {
      // P0 representation
      using P0FEM =
        Dune::PDELab::P0LocalFiniteElementMap<typename Traits::ctype,
                                              typename Traits::Range,
                                              Traits::dim>;
      auto fem =
        std::make_shared<P0FEM>(Dune::GeometryTypes::simplex(Traits::dim));
      using VGFS =
        Dune::PDELab::VectorGridFunctionSpace<typename Traits::EntitySet,
                                              P0FEM,
                                              Traits::dim,
                                              typename Traits::VectorBackend,
                                              typename Traits::VectorBackend,
                                              Dune::PDELab::NoConstraints>;
      auto vgfs = std::make_shared<VGFS>(es, fem);
      vgfs->setDataSetType(
        Dune::PDELab::GridFunctionOutputParameters::Output::cellData);
      using StressVector =
        Dune::PDELab::Backend::Vector<VGFS, typename Traits::ctype>;
      auto stress_container = std::make_shared<StressVector>(vgfs);

      Dune::PDELab::interpolate(stress, *vgfs, *stress_container);
      std::dynamic_pointer_cast<Dune::BlockLab::VisualizationBlock<P, V>>(
        this->parent)
        ->add_dataset(stress_container, "stress_ev_" + std::to_string(_index));
    }
  }

  static std::vector<std::string> blockData()
  {
    auto data = Dune::BlockLab::BlockBase<P, V, i>::blockData();
    // clang-format off
    data.push_back("title: Stress Tensor Eigenvector Visualization   \n"
                   "category: visualization                          \n"
                   "schema:                                          \n"
                   "  index:                                         \n"
                   "    type: integer                                \n"
                   "    default: "
                   + std::to_string(Traits::dim - 1)
                   + "                                  \n"
                   "    min: 0                                       \n"
                   "    max: "
                   + std::to_string(Traits::dim - 1)
                   + "                                       \n"
                   "    meta:                                        \n"
                   "      title: Index of the eigenvector to print   \n"
                   "  continuous:                                    \n"
                   "    type: boolean                                \n"
                   "    default: true                                \n"
                   "    meta:                                        \n"
                   "      title: P1 representation                   \n");
    // clang-format on
    return data;
  }

private:
  size_t _index;
  bool _continuous;
  Material material;
};

#endif
