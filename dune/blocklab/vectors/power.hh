#ifndef DUNE_BLOCKLAB_VECTORS_POWER_HH
#define DUNE_BLOCKLAB_VECTORS_POWER_HH

#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>

#include<memory>
#include<string>


namespace Dune::BlockLab {

  template<typename VectorProvider, std::size_t n>
  class PowerVectorProvider
  {
    using LeafGridFunctionSpace = typename VectorProvider::GridFunctionSpace;
    using VectorBackend = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;

    public:
    // The types of the objects constructed by this provider
    using GridFunctionSpace = Dune::PDELab::PowerGridFunctionSpace<LeafGridFunctionSpace, n, VectorBackend>;
    using Vector = Dune::PDELab::Backend::Vector<GridFunctionSpace, double>;
    using ConstraintsContainer = typename GridFunctionSpace::template ConstraintsContainer<double>::Type;

    // The Parameter types exported by this provider
    using Parameter = typename VectorProvider::Parameter;

    PowerVectorProvider(std::string name,
			std::shared_ptr<VectorProvider> vectorprovider)
      : leafvectorprovider(vectorprovider)
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
        gfs = std::make_shared<GridFunctionSpace>(*leafvectorprovider->getGridFunctionSpace());
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
    std::shared_ptr<VectorProvider> leafvectorprovider;
    std::string name;

    // Internally store the return objects to ensure that the getter methods
    // always return the same object, even if called multiple times
    std::shared_ptr<GridFunctionSpace> gfs;
    std::shared_ptr<Vector> v;
    std::shared_ptr<ConstraintsContainer> cc;
  };

  template<typename std::size_t n, typename VectorProvider>
  auto powerProvider(std::string name, std::shared_ptr<VectorProvider> vectorprovider)
  {
    return std::make_shared<PowerVectorProvider<VectorProvider, n>>(name, vectorprovider);
  }

} // namespace Dune::BlockLab

#endif
