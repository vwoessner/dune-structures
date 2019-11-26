#ifndef DUNE_STRUCTURES_SOLVERSTEPS_INTERPOLATION_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_INTERPOLATION_HH

#include<dune/pdelab.hh>
#include<dune/pdelab/gridfunctionspace/tags.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/typetree/simpletransformationdescriptors.hh>
#include<dune/typetree/utility.hh>

#include<functional>


template<typename SourceNode, typename Transformation>
struct InterpolationLeafNodeTransformation
{
  static const bool recursive = false;

  typedef decltype(Dune::PDELab::makeGridFunctionFromCallable(std::declval<SourceNode>().gridView(),
                                                              std::declval<typename Transformation::Function>()
                                                              )
                   ) transformed_type;

  typedef std::shared_ptr<transformed_type> transformed_storage_type;

  static transformed_type transform(const SourceNode& s, Transformation& t)
  {
    return Dune::PDELab::makeGridFunctionFromCallable(s.gridView(), t.lambdas[t.offset++]);
  }

  static transformed_storage_type transform_storage(std::shared_ptr<const SourceNode> s, Transformation& t)
  {
    return std::make_shared<transformed_type>(Dune::PDELab::makeGridFunctionFromCallable(s->gridView(), t.lambdas[t.offset++]));
  }
};


template<typename RootGFS>
struct GFStoGFTransformation
{
  using Function = std::function<double(typename RootGFS::Traits::GridViewType::template Codim<0>::Entity::Geometry::GlobalCoordinate)>;
  using FunctionArray = std::array<Function, Dune::TypeTree::TreeInfo<RootGFS>::leafCount>;

  GFStoGFTransformation(FunctionArray lambdas)
    : offset(0), lambdas(lambdas)
  {}

  int offset;
  FunctionArray lambdas;
};


template<typename GFS, typename RootGFS>
InterpolationLeafNodeTransformation<GFS, GFStoGFTransformation<RootGFS>>
registerNodeTransformation(GFS* gfs, GFStoGFTransformation<RootGFS>* t, Dune::PDELab::LeafGridFunctionSpaceTag* tag);

template<typename GFS, typename RootGFS>
Dune::TypeTree::SimplePowerNodeTransformation<GFS, GFStoGFTransformation<RootGFS>, Dune::PDELab::PowerGridFunction>
registerNodeTransformation(GFS* gfs, GFStoGFTransformation<RootGFS>* t, Dune::PDELab::PowerGridFunctionSpaceTag* tag);

template<typename GFS, typename RootGFS>
Dune::TypeTree::SimpleCompositeNodeTransformation<GFS, GFStoGFTransformation<RootGFS>, Dune::PDELab::CompositeGridFunction>
registerNodeTransformation(GFS* gfs, GFStoGFTransformation<RootGFS>* t, Dune::PDELab::CompositeGridFunctionSpaceTag* tag);


template<typename Vector>
class InterpolationTransitionStep
  : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;
  using FunctionSignature = typename Base::Range(typename Base::GlobalCoordinate);

  InterpolationTransitionStep(std::function<FunctionSignature> func)
  {
    funcs.fill(func);
  }

  template<typename... FUNCS,
           typename std::enable_if<Dune::TypeTree::TreeInfo<typename Base::GridFunctionSpace>::leafCount == sizeof...(FUNCS), int>::type = 0>
  InterpolationTransitionStep(FUNCS... funcs) : funcs{funcs...}
  {}

  virtual ~InterpolationTransitionStep() {}

  virtual void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer>) override
  {
    auto& gfs = vector->gridFunctionSpace();
    using Trafo = GFStoGFTransformation<typename Base::GridFunctionSpace>;
    Trafo trafo(funcs);

    auto gf = Dune::TypeTree::TransformTree<typename Base::GridFunctionSpace, Trafo>::transform(gfs, trafo);

    std::cout << "Interpolating into solution vector" << std::endl;
    Dune::PDELab::interpolate(gf, gfs, *vector);
  }

  private:
  // Store the lambdas
  std::array<std::function<FunctionSignature>,
             Dune::TypeTree::TreeInfo<typename Base::GridFunctionSpace>::leafCount
             > funcs;
};


#endif
