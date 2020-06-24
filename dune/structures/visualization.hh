#ifndef DUNE_STRUCTURES_VISUALIZATION_HH
#define DUNE_STRUCTURES_VISUALIZATION_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/blocks/enableif.hh>
#include<dune/blocklab/blocks/visualization.hh>
#include<dune/structures/material.hh>
#include<dune/structures/vonmises.hh>

#include<filesystem>
#include<memory>
#include<vector>
#include<variant>


template<typename P, typename V, std::size_t i>
class PhysicalEntityVisualizationBlock
  : public Dune::BlockLab::BlockBase<P, V, i>
{
  public:
  using Traits = Dune::BlockLab::BlockTraits<P, V, i>;

  template<typename Context>
  PhysicalEntityVisualizationBlock(Context&, const Dune::ParameterTree&)
  {}

  virtual ~PhysicalEntityVisualizationBlock() = default;

  struct PhysicalEntityFunction
    : public Dune::VTKFunction<typename Traits::GridView>
  {
    using Base = Dune::VTKFunction<typename Traits::GridView>;
    using Entity = typename Base::Entity;

    PhysicalEntityFunction(typename Traits::GridView gv, std::shared_ptr<std::vector<int>> physical)
      : is(gv.indexSet())
      , physical(physical)
    {}

    virtual ~PhysicalEntityFunction() = default;

    virtual int ncomps() const override
    {
      return 1;
    }

    virtual double evaluate (int comp, const typename Traits::Entity& e,
                             const typename Traits::LocalCoordinate& xi) const override
    {
      return static_cast<double>((*physical)[is.index(e)]);
    }

    virtual std::string name () const override
    {
      return "GMSH Physical Entity Information";
    }

    const typename Traits::GridView::IndexSet& is;
    std::shared_ptr<std::vector<int>> physical;
  };

  virtual void setup() override
  {
    auto gv = this->solver->template getVector<i>()->gridFunctionSpace().gridView();
    auto physical = this->solver->template param<std::shared_ptr<std::vector<int>>>("physical");
    auto function = std::make_shared<PhysicalEntityFunction>(gv, physical);
    std::dynamic_pointer_cast<Dune::BlockLab::VisualizationBlock<P, V, i>>(this->parent)->add_celldata(function);
  }
};


template<typename P, typename V, std::size_t i, typename enabled = Dune::BlockLab::disabled>
class VonMisesStressVisualizationBlock
  : public Dune::BlockLab::DisabledBlock<P, V, i>
{
  public:
  template<typename Context>
  VonMisesStressVisualizationBlock(Context& ctx, const Dune::ParameterTree& config)
   : Dune::BlockLab::DisabledBlock<P, V, i>(ctx, config)
  {}
};


template<typename P, typename V, std::size_t i>
class VonMisesStressVisualizationBlock<P, V, i, Dune::BlockLab::enableBlock<Dune::BlockLab::isDimPower<P, V, i>()>>
  : public Dune::BlockLab::BlockBase<P, V, i>
{
  public:
  using Traits = Dune::BlockLab::BlockTraits<P, V, i>;
  using Material = std::shared_ptr<ElasticMaterialBase<typename Traits::EntitySet, double>>;

  template<typename Context>
  VonMisesStressVisualizationBlock(Context &, const Dune::ParameterTree& config)
    : continuous(config.get<bool>("continuous", false))
  {}

  virtual ~VonMisesStressVisualizationBlock() = default;

  virtual void update_parameter(std::string name, typename Traits::Parameter param) override
  {
    if (name == "material")
      material = std::get<Material>(param);
  }

  virtual void apply() override
  {
    auto vector = this->solver->template getVector<i>();
    auto es = vector->gridFunctionSpace().entitySet();
    auto gfs = vector->gridFunctionSpaceStorage();

    // A grid function for the stress
    VonMisesStressGridFunction<typename Traits::Vector, Traits::dim> stress(*vector, material);

    // Interpolate the stress into a grid function
    if (continuous)
    {
      using P1FEM = Dune::PDELab::PkLocalFiniteElementMap<typename Traits::EntitySet, double, typename Traits::Range, 1>;
      auto p1fem = std::make_shared<P1FEM>(es);
      using P1GFS = Dune::PDELab::GridFunctionSpace<typename Traits::EntitySet, P1FEM, Dune::PDELab::NoConstraints, typename Traits::VectorBackend>;
      auto p1gfs = std::make_shared<P1GFS>(es, p1fem);
      using StressVector = Dune::PDELab::Backend::Vector<P1GFS, typename Traits::ctype>;
      auto stress_container = std::make_shared<StressVector>(p1gfs);

      Dune::PDELab::interpolate(stress, *p1gfs, *stress_container);
      std::dynamic_pointer_cast<Dune::BlockLab::VisualizationBlock<P, V, i>>(this->parent)->add_dataset(stress_container, "vonmises");
    }
    else {
      using P0FEM = Dune::PDELab::P0LocalFiniteElementMap<typename Traits::ctype, typename Traits::Range, Traits::dim>;
      auto p0fem = std::make_shared<P0FEM>(Dune::GeometryTypes::simplex(Traits::dim));
      using P0GFS = Dune::PDELab::GridFunctionSpace<typename Traits::EntitySet, P0FEM, Dune::PDELab::NoConstraints, typename Traits::VectorBackend>;
      auto p0gfs = std::make_shared<P0GFS>(es, p0fem);
      p0gfs->setDataSetType(Dune::PDELab::GridFunctionOutputParameters::Output::cellData);
      using StressVector = Dune::PDELab::Backend::Vector<P0GFS, typename Traits::ctype>;
      auto stress_container = std::make_shared<StressVector>(p0gfs);

      Dune::PDELab::interpolate(stress, *p0gfs, *stress_container);
      std::dynamic_pointer_cast<Dune::BlockLab::VisualizationBlock<P, V, i>>(this->parent)->add_dataset(stress_container, "vonmises");
    }
  }

  private:
  Material material;
  bool continuous;
};


template<typename P, typename V, std::size_t i>
class FibreDistanceVisualizationBlock
  : public Dune::BlockLab::BlockBase<P, V, i>
{
  public:
  using Traits = Dune::BlockLab::BlockTraits<P, V, i>;

  template<typename Context>
  FibreDistanceVisualizationBlock(Context& ctx, const Dune::ParameterTree& config)
    : prestress(ctx.getRootConfig().sub(config.get<std::string>("key")), ctx.getRootConfig())
  {}

  virtual ~FibreDistanceVisualizationBlock() = default;

  virtual void setup() override
  {
    if constexpr (Traits::dim == 3)
    {
      auto vector = this->solver->template getVector<i>();
      auto es = vector->gridFunctionSpace().entitySet();
      using FEM = Dune::PDELab::PkLocalFiniteElementMap<typename Traits::EntitySet, double, typename Traits::Range, 1>;
      auto fem = std::make_shared<FEM>(es);
      using GFS = Dune::PDELab::GridFunctionSpace<typename Traits::EntitySet, FEM, Dune::PDELab::NoConstraints, typename Traits::VectorBackend>;
      auto gfs = std::make_shared<GFS>(es, fem);
      gfs->name("fibredistance");
      using DistanceVector = Dune::PDELab::Backend::Vector<GFS, typename Traits::ctype>;
      auto container = std::make_shared<DistanceVector>(gfs);

      auto lambda = [this](const auto& e, const auto& x){ return this->prestress.distance_to_minimum(e, x); };
      auto gf = Dune::PDELab::makeGridFunctionFromCallable(es.gridView(), lambda);
      Dune::PDELab::interpolate(gf, *gfs, *container);
      std::dynamic_pointer_cast<Dune::BlockLab::VisualizationBlock<P, V, i>>(this->parent)->add_dataset(container, "fibredistance");
    }
  }

  private:
  CurvedFibrePrestress<typename Traits::GridView, double> prestress;
};

#endif
