#ifndef DUNE_BLOCKLAB_VECTORS_PKFEM_HH
#define DUNE_BLOCKLAB_VECTORS_PKFEM_HH

/** A standard PkFEM-based DoF vector provider
 */

#include<dune/pdelab/backend/interface.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/common/partitionviewentityset.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>

#include<memory>
#include<string>


namespace Dune::BlockLab {

  template<typename GridProvider, unsigned int degree>
  class PkFemVectorProvider
  {
    using GridView = typename GridProvider::Grid::LeafGridView;
    using EntitySet = Dune::PDELab::OverlappingEntitySet<GridView>;
    using FiniteElementMap = Dune::PDELab::PkLocalFiniteElementMap<EntitySet, double, double, degree>;
    using ConstraintsAssembler = Dune::PDELab::OverlappingConformingDirichletConstraints;

    public:
    // The types of the objects constructed by this provider
    using GridFunctionSpace = Dune::PDELab::GridFunctionSpace<EntitySet, FiniteElementMap, ConstraintsAssembler>;
    using Vector = Dune::PDELab::Backend::Vector<GridFunctionSpace, double>;
    using ConstraintsContainer = typename GridFunctionSpace::template ConstraintsContainer<double>::Type;

    // The Parameter types exported by this provider
    using Parameter = typename GridProvider::Parameter;

    PkFemVectorProvider(std::shared_ptr<GridProvider> gridprovider,
			std::string name = "solution")
      : gridprovider(gridprovider)
      , name(name)
      , gfs(nullptr)
      , v(nullptr)
      , cc(nullptr)
    {}

    std::string getName()
    {
      return name;
    }

    std::shared_ptr<GridFunctionSpace> getGridFunctionSpace()
    {
      if(!gfs)
      {
        EntitySet es(gridprovider->createGrid()->leafGridView());
        auto fem = std::make_shared<FiniteElementMap>(es);
        gfs = std::make_shared<GridFunctionSpace>(es, fem);
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
    std::shared_ptr<GridProvider> gridprovider;
    std::string name;

    // Internally store the return objects to ensure that the getter methods
    // always return the same object, even if called multiple times
    std::shared_ptr<GridFunctionSpace> gfs;
    std::shared_ptr<Vector> v;
    std::shared_ptr<ConstraintsContainer> cc;
  };

} // namespace Dune::BlockLab

#endif
