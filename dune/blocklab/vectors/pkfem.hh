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


namespace Dune::BlockLab {

  template<typename GridProvider, unsigned int degree>
  class PkFemVectorProvider
  {
    using GridView = typename GridProvider::Grid::LeafGridView;
    using EntitySet = Dune::PDELab::OverlappingEntitySet<GridView>;
    using FiniteElementMap = Dune::PDELab::PkLocalFiniteElementMap<EntitySet, double, double, degree>;
    using ConstraintsAssembler = Dune::PDELab::OverlappingConformingDirichletConstraints;

    public:
    using GridFunctionSpace = Dune::PDELab::GridFunctionSpace<EntitySet, FiniteElementMap, ConstraintsAssembler>;
    using Vector = Dune::PDELab::Backend::Vector<GridFunctionSpace, double>;
    using ConstraintsContainer = typename GridFunctionSpace::template ConstraintsContainer<double>::Type;

    PkFemVectorProvider(std::shared_ptr<GridProvider> gridprovider)
      : gridprovider(gridprovider)
    {}

    std::shared_ptr<GridFunctionSpace> getGridFunctionSpace()
    {
      EntitySet es(gridprovider->createGrid()->leafGridView());
      auto fem = std::make_shared<FiniteElementMap>(es);
      return std::make_shared<GridFunctionSpace>(es, fem);
    }

    std::shared_ptr<Vector> getVector()
    {
      return std::make_shared<Vector>(getGridFunctionSpace());
    }

    std::shared_ptr<ConstraintsContainer> getConstraintsContainer()
    {
      return std::make_shared<ConstraintsContainer>();
    }

    private:
    std::shared_ptr<GridProvider> gridprovider;
  };

} // namespace Dune::BlockLab

#endif
