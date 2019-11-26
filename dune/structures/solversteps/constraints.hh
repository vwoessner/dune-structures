#ifndef DUNE_STRUCTURES_SOLVERSTEPS_CONSTRAINTS_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_CONSTRAINTS_HH


#include<dune/pdelab.hh>
#include<dune/structures/solversteps/base.hh>

#include<functional>


template<typename SourceNode, typename Transformation>
struct ConstraintsLeafNodeTransformation
{
  static const bool recursive = false;

  typedef decltype(Dune::PDELab::makeBoundaryConditionFromCallable(std::declval<SourceNode>().gridView(),
                                                                   std::declval<typename Transformation::Function>()
                                                                   )
                   ) transformed_type;

  typedef std::shared_ptr<transformed_type> transformed_storage_type;

  static transformed_type transform(const SourceNode& s, Transformation& t)
  {
    return Dune::PDELab::makeBoundaryConditionFromCallable(s.gridView(), t.lambdas[t.offset++]);
  }

  static transformed_storage_type transform_storage(std::shared_ptr<const SourceNode> s, Transformation& t)
  {
    return std::make_shared<transformed_type>(Dune::PDELab::makeBoundaryConditionFromCallable(s->gridView(), t.lambdas[t.offset++]));
  }
};


template<typename RootGFS>
struct GFStoConstraintsTransformation
{
  using Function = std::function<bool(typename RootGFS::Traits::GridViewType::template Codim<0>::Entity::Geometry::GlobalCoordinate)>;
  using FunctionArray = std::array<Function, Dune::TypeTree::TreeInfo<RootGFS>::leafCount>;

  GFStoConstraintsTransformation(FunctionArray lambdas)
    : offset(0), lambdas(lambdas)
  {}

  int offset;
  FunctionArray lambdas;
};


template<typename GFS, typename RootGFS>
ConstraintsLeafNodeTransformation<GFS, GFStoConstraintsTransformation<RootGFS>>
registerNodeTransformation(GFS* gfs, GFStoConstraintsTransformation<RootGFS>* t, Dune::PDELab::LeafGridFunctionSpaceTag* tag);

template<typename GFS, typename RootGFS>
Dune::TypeTree::SimplePowerNodeTransformation<GFS, GFStoConstraintsTransformation<RootGFS>, Dune::PDELab::PowerConstraintsParameters>
registerNodeTransformation(GFS* gfs, GFStoConstraintsTransformation<RootGFS>* t, Dune::PDELab::PowerGridFunctionSpaceTag* tag);

template<typename GFS, typename RootGFS>
Dune::TypeTree::SimpleCompositeNodeTransformation<GFS, GFStoConstraintsTransformation<RootGFS>, Dune::PDELab::CompositeConstraintsParameters>
registerNodeTransformation(GFS* gfs, GFStoConstraintsTransformation<RootGFS>* t, Dune::PDELab::CompositeGridFunctionSpaceTag* tag);


template<typename Vector>
class ConstraintsTransitionStep : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;
  using FunctionSignature = bool(typename Base::GlobalCoordinate);

  ConstraintsTransitionStep(std::function<FunctionSignature> func)
  {
    funcs.fill(func);
  }

  template<typename... FUNCS,
           typename std::enable_if<Dune::TypeTree::TreeInfo<typename Base::GridFunctionSpace>::leafCount == sizeof...(FUNCS), int>::type = 0>
  ConstraintsTransitionStep(FUNCS... funcs) : funcs{funcs...}
  {}

  virtual ~ConstraintsTransitionStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> constraintscontainer) override
  {
    auto& gfs = vector->gridFunctionSpace();
    using Trafo = GFStoConstraintsTransformation<typename Base::GridFunctionSpace>;
    Trafo trafo(funcs);

    auto bctype = Dune::TypeTree::TransformTree<typename Base::GridFunctionSpace, Trafo>::transform(gfs, trafo);
    std::cout << "Assembling constraints!" << std::endl;
    Dune::PDELab::constraints(bctype, gfs, *constraintscontainer);
  }

  private:
  // Store the lambdas
  std::array<std::function<FunctionSignature>,
             Dune::TypeTree::TreeInfo<typename Base::GridFunctionSpace>::leafCount
             > funcs;
};

#endif
