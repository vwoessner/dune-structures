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


template<typename Vector>
class VisualizationStepBase
  : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  VisualizationStepBase()
    : vtkwriter(0)
  {}

  virtual ~VisualizationStepBase() {}

  virtual void set_vtkwriter(std::shared_ptr<Dune::VTKWriter<typename Base::GridView>> vtkwriter_)
  {
    vtkwriter = vtkwriter_;
  }

  protected:
  std::shared_ptr<Dune::VTKWriter<typename Base::GridView>> vtkwriter;
};


template<typename Vector>
class VisualizationStep
  : public VisualizationStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  VisualizationStep()
    : steps(0)
  {}

  virtual ~VisualizationStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    // Instantiate a VTKWriter instance
    auto gv = vector->gridFunctionSpace().gridView();
    auto writer = std::make_shared<Dune::VTKWriter<typename Base::GridView>>(gv);

    // Make it known to all child steps
    this->set_vtkwriter(writer);

    for (auto step: steps)
      step->apply(vector, cc);

    this->vtkwriter->write("output", Dune::VTK::ascii);
  }

  virtual void set_vtkwriter(std::shared_ptr<Dune::VTKWriter<typename Base::GridView>> vtkwriter_) override
  {
    // Set it in this instance
    this->vtkwriter = vtkwriter_;

    for (auto step : steps)
      step->set_vtkwriter(vtkwriter_);
  }

  void add(std::shared_ptr<VisualizationStepBase<Vector>> step)
  {
    steps.push_back(step);
  }

  template<typename STEP>
  void add(STEP& step)
  {
    add(Dune::stackobject_to_shared_ptr(step));
  }

  private:
  std::vector<std::shared_ptr<VisualizationStepBase<Vector>>> steps;
};


template<typename Vector>
class SolutionVisualizationStep
  : public VisualizationStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  virtual ~SolutionVisualizationStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    const auto& gfs = vector->gridFunctionSpace();
    auto ugly_gfs_pointer = Dune::stackobject_to_shared_ptr(gfs);
    Dune::PDELab::addSolutionToVTKWriter(*(this->vtkwriter), ugly_gfs_pointer, vector);
  }
};


template<typename Vector>
class MPIRankVisualizationStep
  : public VisualizationStepBase<Vector>
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

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    RankDummyContainer container(helper, vector->gridFunctionSpace().gridView());
    this->vtkwriter->addCellData(container, "mpirank");
  }

  private:
  Dune::MPIHelper& helper;
};


template<typename Vector>
class PhysicalEntityVisualizationStep
  : public VisualizationStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  PhysicalEntityVisualizationStep(std::shared_ptr<std::vector<int>> physical)
    : physical(physical)
  {}

  virtual ~PhysicalEntityVisualizationStep() {}

  virtual void apply(std::shared_ptr<Vector>, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    this->vtkwriter->addCellData(*physical, "gmshPhysical");
  }

  private:
  std::shared_ptr<std::vector<int>> physical;
};


template<typename Vector>
class VonMisesStressVisualizationStep
  : public VisualizationStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  VonMisesStressVisualizationStep(std::shared_ptr<ElasticMaterialBase<typename Base::EntitySet, typename Base::Range>> material)
    : material(material)
  {}

  virtual ~VonMisesStressVisualizationStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    auto es = vector->gridFunctionSpace().entitySet();

    // A grid function for the stress
    VonMisesStressGridFunction stress(*vector, material);

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
  std::shared_ptr<ElasticMaterialBase<typename Base::EntitySet, typename Base::Range>> material;
};

#endif
