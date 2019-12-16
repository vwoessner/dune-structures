#ifndef DUNE_STRUCTURES_SOLVERSTEPS_VISUALIZATION_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_VISUALIZATION_HH

#include<dune/common/shared_ptr.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/pdelab.hh>
#include<dune/structures/material.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/vonmises.hh>

#include<memory>
#include<vector>


template<typename Vector, bool instationary>
struct VTKWriterChooser
{
  using GV = typename Vector::GridFunctionSpace::Traits::GridViewType;
  using type = Dune::VTKWriter<GV>;
};

template<typename Vector>
struct VTKWriterChooser<Vector, true>
{
  using GV = typename Vector::GridFunctionSpace::Traits::GridViewType;
  using type = Dune::VTKSequenceWriter<GV>;
};




template<typename Vector, bool instationary>
class VisualizationStepBase
  : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;
  using VTKWriter = typename VTKWriterChooser<Vector, instationary>::type;

  VisualizationStepBase()
    : vtkwriter(0)
  {}

  virtual ~VisualizationStepBase() {}

  virtual void set_vtkwriter(std::shared_ptr<VTKWriter> vtkwriter_)
  {
    vtkwriter = vtkwriter_;
  }

  protected:
  std::shared_ptr<VTKWriter> vtkwriter;
};


template<typename Vector, bool instationary=false>
class VisualizationStep
  : public VisualizationStepBase<Vector, instationary>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;
  using VTKWriter = typename VTKWriterChooser<Vector, instationary>::type;

  VisualizationStep(std::string name="output")
    : steps(0), time(0.0), name(name), path(""), extendpath("")
  {}

  VisualizationStep(const std::string& name,
                    const std::string& path,
                    const std::string& extendpath = "")
    : steps(0), time(0.0), name(name), path(path), extendpath(extendpath)
  {}

  VisualizationStep(const Dune::ParameterTree& config)
    : steps(0), time(0.0)
    , name(config.get<std::string>("name"))
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

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    // Instantiate a VTKWriter instance
    auto gv = vector->gridFunctionSpace().gridView();

    if constexpr (instationary)
      this->set_vtkwriter(std::make_shared<typename VTKWriterChooser<Vector, true>::type>(std::make_shared<typename VTKWriterChooser<Vector, false>::type>(gv), name, path, extendpath));
    else
      this->set_vtkwriter(std::make_shared<typename VTKWriterChooser<Vector, false>::type>(gv));

    for (auto step: steps)
      step->pre(vector, cc);
  }

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    for (auto step: steps)
      step->apply(vector, cc);

    if constexpr (instationary)
      this->vtkwriter->write(time, Dune::VTK::appendedraw);
    else
      this->vtkwriter->write("output", Dune::VTK::ascii);
  }

  virtual void set_vtkwriter(std::shared_ptr<VTKWriter> vtkwriter_) override
  {
    // Set it in this instance
    this->vtkwriter = vtkwriter_;

    for (auto step : steps)
      step->set_vtkwriter(vtkwriter_);
  }

  void add(std::shared_ptr<VisualizationStepBase<Vector, instationary>> step)
  {
    steps.push_back(step);
  }

  template<typename STEP>
  void add(STEP& step)
  {
    add(Dune::stackobject_to_shared_ptr(step));
  }

  private:
  std::vector<std::shared_ptr<VisualizationStepBase<Vector, instationary>>> steps;
  double time;
  std::string name;
  std::string path;
  std::string extendpath;
};


template<typename Vector, bool instationary=false>
class SolutionVisualizationStep
  : public VisualizationStepBase<Vector, instationary>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  virtual ~SolutionVisualizationStep() {}

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    if constexpr (instationary)
    {
      auto gfs = vector->gridFunctionSpaceStorage();
      auto sgfs1 = std::make_shared<Dune::PDELab::GridFunctionSubSpace<typename Base::GridFunctionSpace, Dune::TypeTree::TreePath<0>>>(gfs);
      Dune::PDELab::addSolutionToVTKWriter(*(this->vtkwriter), sgfs1, vector);
      auto sgfs2 = std::make_shared<Dune::PDELab::GridFunctionSubSpace<typename Base::GridFunctionSpace, Dune::TypeTree::TreePath<1>>>(gfs);
      Dune::PDELab::addSolutionToVTKWriter(*(this->vtkwriter), sgfs2, vector);
    }
    else
    {
      auto gfs = vector->gridFunctionSpaceStorage();
      Dune::PDELab::addSolutionToVTKWriter(*(this->vtkwriter), gfs, vector);
    }
  }
};


template<typename Vector, bool instationary=false>
class MPIRankVisualizationStep
  : public VisualizationStepBase<Vector, instationary>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

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

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    RankDummyContainer container(helper, vector->gridFunctionSpace().gridView());
    this->vtkwriter->addCellData(container, "mpirank");
  }

  private:
  Dune::MPIHelper& helper;
};


template<typename Vector, bool instationary=false>
class PhysicalEntityVisualizationStep
  : public VisualizationStepBase<Vector, instationary>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  PhysicalEntityVisualizationStep(std::shared_ptr<std::vector<int>> physical)
    : physical(physical)
  {}

  virtual ~PhysicalEntityVisualizationStep() {}

  virtual void pre(std::shared_ptr<Vector>, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    this->vtkwriter->addCellData(*physical, "gmshPhysical");
  }

  private:
  std::shared_ptr<std::vector<int>> physical;
};


template<typename Vector, bool instationary=false>
class VonMisesStressVisualizationStep
  : public VisualizationStepBase<Vector, instationary>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  virtual ~VonMisesStressVisualizationStep() {}

  virtual void update_parameter(std::string name, typename Base::Parameter param) override
  {
    if (name == "material")
      material = std::get<std::shared_ptr<typename Base::Material>>(param);
  }

  virtual void pre(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    auto es = vector->gridFunctionSpace().entitySet();

    // A grid function for the stress
    VonMisesStressGridFunction<Vector> stress(*vector, material);

    // Interpolate the stress into a grid function
    using P0FEM = Dune::PDELab::P0LocalFiniteElementMap<typename Base::ctype, typename Base::Range, 3>;
    auto p0fem = std::make_shared<P0FEM>(Dune::GeometryTypes::simplex(3));
    using P0GFS = Dune::PDELab::GridFunctionSpace<typename Base::EntitySet, P0FEM, Dune::PDELab::NoConstraints, typename Base::VectorBackend>;
    auto p0gfs = std::make_shared<P0GFS>(es, p0fem);
    p0gfs->name("vonmises");
    p0gfs->setDataSetType(Dune::PDELab::GridFunctionOutputParameters::Output::cellData);
    using StressVector = Dune::PDELab::Backend::Vector<P0GFS, typename Base::ctype>;
    auto stress_container = std::make_shared<StressVector>(p0gfs);

    Dune::PDELab::interpolate(stress, *p0gfs, *stress_container);
    Dune::PDELab::addSolutionToVTKWriter(*(this->vtkwriter), p0gfs, stress_container);
  }

  private:
  std::shared_ptr<typename Base::Material> material;
};

#endif
