#ifndef DUNE_BLOCKLAB_VECTORS_COMPOSITE_HH
#define DUNE_BLOCKLAB_VECTORS_COMPOSITE_HH

#include<dune/blocklab/utilities/tuplecat.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/gridfunctionspace/compositegridfunctionspace.hh>

#include<memory>
#include<string>
#include<tuple>


namespace Dune::BlockLab {

  template<typename... VectorProvider>
  class CompositeVectorProvider
  {
    using VectorBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
    using OrderingTag = Dune::PDELab::LexicographicOrderingTag;

    public:
    // The types of the objects constructed by this provider
    using GridFunctionSpace = Dune::PDELab::CompositeGridFunctionSpace<VectorBackend, OrderingTag, typename VectorProvider::GridFunctionSpace...>;
    using Vector = Dune::PDELab::Backend::Vector<GridFunctionSpace, double>;
    using ConstraintsContainer = typename GridFunctionSpace::template ConstraintsContainer<double>::Type;

    // The Parameter types exported by this provider
    using Parameter = tuple_cat_t<typename VectorProvider::Parameter...>;

    CompositeVectorProvider(std::string name,
			    std::shared_ptr<VectorProvider>... vectorprovider)
      : leafvectorprovider(vectorprovider...)
      , name(name)
      , gfs(nullptr)
      , v(nullptr)
      , cc(nullptr)
    {}

    std::string getName()
    {
      return name;
    }

    std::shared_ptr<const GridFunctionSpace> getGridFunctionSpace()
    {
      if(!gfs)
      {
        gfs = std::apply([](auto... v) {
                           return std::make_shared<GridFunctionSpace>(v->getGridFunctionSpace()...);
                         },
			 leafvectorprovider);
        gfs->name(name);
      }
      return gfs;
    }

    std::shared_ptr<Vector> getVector()
    {
      if (!v)
	v = std::make_shared<Vector>(getGridFunctionSpace());
      return v;
    }

    std::shared_ptr<ConstraintsContainer> getConstraintsContainer()
    {
      if (!cc)
        cc = std::make_shared<ConstraintsContainer>();
      return cc;
    }

    private:
    std::tuple<std::shared_ptr<VectorProvider>...> leafvectorprovider;
    std::string name;

    // Internally store the return objects to ensure that the getter methods
    // always return the same object, even if called multiple times
    std::shared_ptr<GridFunctionSpace> gfs;
    std::shared_ptr<Vector> v;
    std::shared_ptr<ConstraintsContainer> cc;
  };

  template<typename... VectorProvider>
  auto compositeProvider(std::string name, std::shared_ptr<VectorProvider>... vectorprovider)
  {
    return std::make_shared<CompositeVectorProvider<VectorProvider...>>(name, vectorprovider...);
  }

} // namespace Dune::BlockLab

#endif
