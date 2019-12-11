#ifndef DUNE_STRUCTURES_CALLABLEADAPTERS_HH
#define DUNE_STRUCTURES_CALLABLEADAPTERS_HH

#include<dune/pdelab.hh>
#include<dune/typetree/simpletransformationdescriptors.hh>
#include<dune/typetree/utility.hh>

#include<memory>


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

  template<typename... FUNCS>
  GFStoGFTransformation(FUNCS... lambdas)
    : offset(0), lambdas{lambdas...}
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


template<typename GFS, typename... FUNCS>
auto makeGridFunctionTreeFromCallables(const GFS& gfs, FUNCS... funcs)
{
  using Trafo = GFStoGFTransformation<GFS>;
  Trafo trafo(funcs...);
  return Dune::TypeTree::TransformTree<GFS, Trafo>::transform(gfs, trafo);
}



template<typename P, typename SourceNode, typename Transformation>
struct InstationaryInterpolationLeafNodeTransformation
{
  static const bool recursive = false;

  typedef decltype(Dune::PDELab::makeInstationaryGridFunctionFromCallable(std::declval<SourceNode>().gridView(),
                                                                          std::declval<typename Transformation::Function>(),
                                                                          std::declval<P&>()
                                                                          )
                   ) transformed_type;

  typedef std::shared_ptr<transformed_type> transformed_storage_type;

  static transformed_type transform(const SourceNode& s, Transformation& t)
  {
    return Dune::PDELab::makeInstationaryGridFunctionFromCallable(s.gridView(), t.lambdas[t.offset++], t.tc);
  }

  static transformed_storage_type transform_storage(std::shared_ptr<const SourceNode> s, Transformation& t)
  {
    return std::make_shared<transformed_type>(Dune::PDELab::makeInstationaryGridFunctionFromCallable(s->gridView(), t.lambdas[t.offset++], t.tc));
  }
};


template<typename P, typename RootGFS>
struct GFStoInstationaryGFTransformation
{
  using Function = std::function<double(typename RootGFS::Traits::GridViewType::template Codim<0>::Entity::Geometry::GlobalCoordinate)>;
  using FunctionArray = std::array<Function, Dune::TypeTree::TreeInfo<RootGFS>::leafCount>;

  template<typename... FUNCS>
  GFStoInstationaryGFTransformation(P& tc, FUNCS... lambdas)
    : tc(tc), offset(0), lambdas{lambdas...}
  {}

  P& tc;
  int offset;
  FunctionArray lambdas;
};


template<typename P, typename GFS, typename RootGFS>
InstationaryInterpolationLeafNodeTransformation<P, GFS, GFStoInstationaryGFTransformation<P, RootGFS>>
registerNodeTransformation(GFS* gfs, GFStoInstationaryGFTransformation<P, RootGFS>* t, Dune::PDELab::LeafGridFunctionSpaceTag* tag);

template<typename P, typename GFS, typename RootGFS>
Dune::TypeTree::SimplePowerNodeTransformation<GFS, GFStoInstationaryGFTransformation<P, RootGFS>, Dune::PDELab::PowerGridFunction>
registerNodeTransformation(GFS* gfs, GFStoInstationaryGFTransformation<P, RootGFS>* t, Dune::PDELab::PowerGridFunctionSpaceTag* tag);

template<typename P, typename GFS, typename RootGFS>
Dune::TypeTree::SimpleCompositeNodeTransformation<GFS, GFStoInstationaryGFTransformation<P, RootGFS>, Dune::PDELab::CompositeGridFunction>
registerNodeTransformation(GFS* gfs, GFStoInstationaryGFTransformation<P, RootGFS>* t, Dune::PDELab::CompositeGridFunctionSpaceTag* tag);


template<typename P, typename GFS, typename... FUNCS>
auto makeInstationaryGridFunctionTreeFromCallables(P& tc,
                                                   const GFS& gfs,
                                                   FUNCS... funcs)
{
  using Trafo = GFStoInstationaryGFTransformation<P, GFS>;
  Trafo trafo(tc, funcs...);
  return Dune::TypeTree::TransformTree<GFS, Trafo>::transform(gfs, trafo);
}


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


template<typename GFS, typename... FUNCS>
auto makeBoundaryConditionTreeFromCallables(const GFS& gfs, FUNCS... funcs)
{
  using Trafo = GFStoConstraintsTransformation<GFS>;
  Trafo trafo(funcs...);
  return Dune::TypeTree::TransformTree<GFS, Trafo>::transform(gfs, trafo);
}

#endif