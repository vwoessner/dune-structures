#ifndef DUNE_STRUCTURES_SOLVERSTEPS_ELASTICITY_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_ELASTICITY_HH

#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/newton.hh>

#include"elasticity_operator.hh"


template<typename Vector>
class ElasticitySolverStep
  : public NewtonSolverTransitionStep<Vector,
                                      ElasticityOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                                         typename TransitionSolverStepBase<Vector>::GridFunctionSpace>
                                      >
{
  public:
  using Base = TransitionSolverStepBase<Vector>;
  using LocalOperator = ElasticityOperator<typename TransitionSolverStepBase<Vector>::GridFunctionSpace,
                                           typename TransitionSolverStepBase<Vector>::GridFunctionSpace>;

  ElasticitySolverStep(std::shared_ptr<std::vector<int>> physical,
                       const Dune::ParameterTree& params
                       )
    : physical(physical),
      params(params),
      material_params(params.sub("material"))
  {}

  virtual ~ElasticitySolverStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    // Parse the material
    auto& gfs = vector->gridFunctionSpace();
    auto es = gfs.entitySet();
    auto material = parse_material<double>(es, physical, material_params);

    // Instantiate the local operator
    this->localoperator = std::make_shared<LocalOperator>(gfs, gfs, params, material);

    // Call the base class
    this->NewtonSolverTransitionStep<Vector, LocalOperator>::apply(vector, cc);
  }

  private:
  std::shared_ptr<std::vector<int>> physical;
  const Dune::ParameterTree& params;
  Dune::ParameterTree material_params;
};


template<typename ES, typename RangeType = double>
auto elasticity_setup(ES es)
{
  // Set up finite element maps...
  using FEM = Dune::PDELab::PkLocalFiniteElementMap<ES, double, RangeType, 1>;
  auto fem = std::make_shared<FEM>(es);

  // Set up grid function spaces...
  using VB = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
  using CASS = Dune::PDELab::ConformingDirichletConstraints;
  using GFS = Dune::PDELab::VectorGridFunctionSpace<ES, FEM, 3, VB, VB, CASS>;
  auto gfs = std::make_shared<GFS>(es, fem);
  gfs->name("displacement");
  gfs->update();
  std::cout << "Set up a grid function space with " << gfs->size() << " dofs!" << std::endl;

  // Setting up constraints container
  using CC = typename GFS::template ConstraintsContainer<RangeType>::Type;
  auto cc = std::make_shared<CC>();
  cc->clear();

  // Setting up container
  using V = Dune::PDELab::Backend::Vector<GFS, double>;
  auto x = std::make_shared<V>(gfs);

  return std::make_pair(x, cc);
}

#endif
