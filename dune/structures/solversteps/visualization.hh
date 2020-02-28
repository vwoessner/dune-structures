#ifndef DUNE_STRUCTURES_SOLVERSTEPS_VISUALIZATION_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_VISUALIZATION_HH

#include<dune/common/shared_ptr.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/pdelab.hh>
#include<dune/structures/material.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/vonmises.hh>

#include<filesystem>
#include<memory>
#include<vector>
#include<variant>


template<typename Vector, bool instationary>
struct VTKWriterChooser
{
  using GV = typename Vector::GridFunctionSpace::Traits::GridViewType;
  using type = Dune::VTKWriter<GV>;
  using ptype = std::shared_ptr<type>;
};


template<typename Vector>
struct VTKWriterChooser<Vector, true>
{
  using GV = typename Vector::GridFunctionSpace::Traits::GridViewType;
  using type = Dune::VTKSequenceWriter<GV>;
  using ptype = std::shared_ptr<type>;
};


template<typename... V>
class VisualizationStep;


template<typename... V>
class VisualizationStepBase
  : public TransitionSolverStepBase<V...>
{
  public:
  using Base = TransitionSolverStepBase<V...>;

  void set_parent(VisualizationStep<V...>* parent_)
  {
    parent = parent_;
  }

  protected:
  VisualizationStep<V...>* parent;
};


template<typename... V>
class VisualizationStep
  : public StepCollectionStep<V...>
{
  public:
  using Base = TransitionSolverStepBase<V...>;
  using VTKWriter = std::variant<
      std::shared_ptr<typename VTKWriterChooser<typename Base::Vector, true>::type>,
      std::shared_ptr<typename VTKWriterChooser<typename Base::Vector, false>::type>>;

  VisualizationStep(bool instationary = false,
                    std::string name="output")
    : instationary(instationary), time(0.0), name(name), path(""), extendpath("")
  {}

  VisualizationStep(bool instationary,
                    const std::string& name,
                    const std::string& path,
                    const std::string& extendpath = "")
    : instationary(instationary)
    , time(0.0)
    , name(name)
    , path(path)
    , extendpath(extendpath)
  {}

  VisualizationStep(const Dune::ParameterTree& config)
    : instationary(config.get<bool>("instationary", false))
    , time(0.0)
    , name(config.get<std::string>("name", "output"))
    , path(config.get<std::string>("path", ""))
    , extendpath(config.get<std::string>("extendpath", ""))
  {}

  virtual ~VisualizationStep() {}

  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    if (name == "time")
      time = std::get<double>(param);

    for (auto step : this->steps)
      step->update_parameter(name, param);
  }

  virtual void pre(std::shared_ptr<typename Base::Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    // Instantiate a VTKWriter instance
    auto gv = vector->gridFunctionSpace().gridView();

    if (instationary)
      vtkwriter = std::make_shared<typename VTKWriterChooser<typename Base::Vector, true>::type>(std::make_shared<typename VTKWriterChooser<typename Base::Vector, false>::type>(gv), name, path, extendpath);
    else
      vtkwriter = std::make_shared<typename VTKWriterChooser<typename Base::Vector, false>::type>(gv);

    for (auto step: this->steps)
    {
      auto vsp = dynamic_cast<VisualizationStepBase<V...>*>(step.get());
      if (vsp)
        vsp->set_parent(this);
      step->pre(vector, cc);
    }
  }

  virtual void apply(std::shared_ptr<typename Base::Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    for (auto step: this->steps)
      step->apply(vector, cc);

    if (instationary)
    {
      std::filesystem::create_directory(std::filesystem::current_path().append(path));
      std::get<typename VTKWriterChooser<typename Base::Vector, true>::ptype>(vtkwriter)->write(time, Dune::VTK::appendedraw);
    }
    else
      std::get<typename VTKWriterChooser<typename Base::Vector, false>::ptype>(vtkwriter)->write(name, Dune::VTK::ascii);
  }

  template<typename Container>
  void add_dataset(std::shared_ptr<Container> container)
  {
    if (instationary)
      Dune::PDELab::addSolutionToVTKWriter(
          *(std::get<typename VTKWriterChooser<typename Base::Vector, true>::ptype>(vtkwriter)),
          container->gridFunctionSpaceStorage(),
          container);
    else
      Dune::PDELab::addSolutionToVTKWriter(
          *(std::get<typename VTKWriterChooser<typename Base::Vector, false>::ptype>(vtkwriter)),
          container->gridFunctionSpaceStorage(),
          container);
  };

  template<typename Container>
  void add_celldata(std::shared_ptr<Container> container, std::string name)
  {
    if (instationary)
      std::get<typename VTKWriterChooser<typename Base::Vector, true>::ptype>(vtkwriter)->addCellData(*container, name);
    else
      std::get<typename VTKWriterChooser<typename Base::Vector, false>::ptype>(vtkwriter)->addCellData(*container, name);
  }

  private:
  bool instationary;
  double time;
  std::string name;
  std::string path;
  std::string extendpath;
  VTKWriter vtkwriter;
};


template<typename... V>
class SolutionVisualizationStep
  : public VisualizationStepBase<V...>
{
  public:
  using Base = TransitionSolverStepBase<V...>;

  virtual ~SolutionVisualizationStep() {}

  virtual void pre(std::shared_ptr<typename Base::Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    this->parent->add_dataset(vector);
  }
};


template<typename... V>
class MPIRankVisualizationStep
  : public VisualizationStepBase<V...>
{
  public:
  using Base = TransitionSolverStepBase<V...>;

  struct RankDummyContainer
  {
    RankDummyContainer(Dune::MPIHelper& helper, typename Base::GridView gv) : rank(helper.rank()), size_(gv.size(0))
    {}

    double operator[](std::size_t i) const
    {
      return rank;
    }

    std::size_t size() const
    {
      return size_;
    }

    double rank;
    std::size_t size_;
  };

  MPIRankVisualizationStep(Dune::MPIHelper& helper)
    : helper(helper)
  {}

  virtual ~MPIRankVisualizationStep() {}

  virtual void pre(std::shared_ptr<typename Base::Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    auto container = std::make_shared<RankDummyContainer>(helper, vector->gridFunctionSpace().gridView());
    this->parent->add_celldata(container, "mpirank");
  }

  private:
  Dune::MPIHelper& helper;
};


template<typename... V>
class PhysicalEntityVisualizationStep
  : public VisualizationStepBase<V...>
{
  public:
  using Base = TransitionSolverStepBase<V...>;

  PhysicalEntityVisualizationStep(std::shared_ptr<std::vector<int>> physical)
    : physical(physical)
  {}

  virtual ~PhysicalEntityVisualizationStep() {}

  virtual void pre(std::shared_ptr<typename Base::Vector>, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    this->parent->add_celldata(physical, "gmshPhysical");
  }

  private:
  std::shared_ptr<std::vector<int>> physical;
};


template<typename... V>
class VonMisesStressVisualizationStep
  : public VisualizationStepBase<V...>
{
  public:
  using Base = TransitionSolverStepBase<V...>;

  virtual ~VonMisesStressVisualizationStep() {}

  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    if (name == "material")
      material = std::get<std::shared_ptr<typename Base::Material>>(param);
  }

  virtual void apply(std::shared_ptr<typename Base::Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    auto es = vector->gridFunctionSpace().entitySet();

    // A grid function for the stress
    VonMisesStressGridFunction<typename Base::Vector> stress(*vector, material);

    // Interpolate the stress into a grid function
    using P0FEM = Dune::PDELab::P0LocalFiniteElementMap<typename Base::ctype, typename Base::Range, Base::dim>;
    auto p0fem = std::make_shared<P0FEM>(Dune::GeometryTypes::simplex(Base::dim));
    using P0GFS = Dune::PDELab::GridFunctionSpace<typename Base::EntitySet, P0FEM, Dune::PDELab::NoConstraints, typename Base::VectorBackend>;
    auto p0gfs = std::make_shared<P0GFS>(es, p0fem);
    p0gfs->name("vonmises");
    p0gfs->setDataSetType(Dune::PDELab::GridFunctionOutputParameters::Output::cellData);
    using StressVector = Dune::PDELab::Backend::Vector<P0GFS, typename Base::ctype>;
    auto stress_container = std::make_shared<StressVector>(p0gfs);

    Dune::PDELab::interpolate(stress, *p0gfs, *stress_container);
    this->parent->add_dataset(stress_container);
  }

  private:
  std::shared_ptr<typename Base::Material> material;
};


template<typename... V>
class FibreDistanceVisualizationStep
  : public VisualizationStepBase<V...>
{
  public:
  using Base = TransitionSolverStepBase<V...>;

  FibreDistanceVisualizationStep(const Dune::ParameterTree& params, const Dune::ParameterTree& rootparams)
    : prestress(rootparams.sub(params.get<std::string>("key")), rootparams)
  {}

  virtual ~FibreDistanceVisualizationStep() {}

  virtual void pre(std::shared_ptr<typename Base::Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    if constexpr (Base::dim == 3)
    {
      auto es = vector->gridFunctionSpace().entitySet();
      using FEM = Dune::PDELab::PkLocalFiniteElementMap<typename Base::EntitySet, double, typename Base::Range, 1>;
      auto fem = std::make_shared<FEM>(es);
      using GFS = Dune::PDELab::GridFunctionSpace<typename Base::EntitySet, FEM, Dune::PDELab::NoConstraints, typename Base::VectorBackend>;
      auto gfs = std::make_shared<GFS>(es, fem);
      gfs->name("fibredistance");
      using DistanceVector = Dune::PDELab::Backend::Vector<GFS, typename Base::ctype>;
      auto container = std::make_shared<DistanceVector>(gfs);

      auto lambda = [this](const auto& e, const auto& x){ return this->prestress.distance_to_minimum(e, x); };
      auto gf = Dune::PDELab::makeGridFunctionFromCallable(es.gridView(), lambda);
      Dune::PDELab::interpolate(gf, *gfs, *container);
      this->parent->add_dataset(container);
    }
  }

  private:
  CurvedFibrePrestress<typename Base::GridView, double> prestress;
};

#endif
