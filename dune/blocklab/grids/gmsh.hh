#ifndef DUNE_BLOCKLAB_GRIDS_GMSH_HH
#define DUNE_BLOCKLAB_GRIDS_GMSH_HH

#include<dune/common/parametertree.hh>
#include<dune/grid/common/gridfactory.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/uggrid.hh>
#include<dune/pdelab/common/partitionviewentityset.hh>

#include<map>
#include<memory>
#include<tuple>
#include<vector>


namespace Dune::BlockLab {

  namespace impl {

    template<typename IDSet, typename DataMap>
    class PhysicalInfoDataHandle
      : public Dune::CommDataHandleIF<PhysicalInfoDataHandle<IDSet, DataMap>, int>
    {
      public:
      PhysicalInfoDataHandle(const IDSet& idset, DataMap& datamap)
        : idset(idset), datamap(datamap)
      {}

      bool contains (int dim, int codim) const
      {
        return codim == 0;
      }

      bool fixedSize (int dim, int codim) const
      {
        return true;
      }

      template<class EntityType>
      std::size_t size (const EntityType& e) const
      {
        return 1;
      }

      template<class MessageBufferImp, class EntityType>
      void gather (MessageBufferImp& buff, const EntityType& e) const
      {
        buff.write(datamap[idset.id(e)]);
      }

      template<class MessageBufferImp, class EntityType>
      void scatter (MessageBufferImp& buff, const EntityType& e, size_t n)
      {
        buff.read(datamap[idset.id(e)]);
      }

      private:
      const IDSet& idset;
      DataMap& datamap;
    };

  } // namespace impl

  /** We only provide a specialization for Dune::UGGrid here. This has the following reason:
   *  The way we parse the data for physical entities and use them later on makes the silent
   *  assumption that the grid factory insertion order and the level 0 index set are the same.
   *  This holds for UG, but it is an implementation detail, not an interface property.
   *  The code should query the factory for the insertion indices and reorder the data
   *  to match the level 0 index. Once that is done, below specialization can become the
   *  implementation.
   */
  template<typename GridImpl>
  class GMSHGridProvider
  {};

  template<int dimension>
  class GMSHGridProvider<Dune::UGGrid<dimension>>
  {
    public:
    using Grid = Dune::UGGrid<dimension>;
    using Parameter = std::tuple<std::shared_ptr<std::vector<int>>>;
    using EntitySet = Dune::PDELab::OverlappingEntitySet<typename Grid::LeafGridView>;
    static constexpr int dim = Grid::dimension;

    GMSHGridProvider(const Dune::ParameterTree& config)
      : config(config)
      , physical(nullptr)
    {}

    std::shared_ptr<Grid> createGrid()
    {
      if(!grid)
      {
        auto mshfile = config.get<std::string>("filename");
	Dune::GridFactory<Grid> factory;
	std::vector<int> boundary;

        // Read the grid (only operative on Rank 0)
        Dune::GmshReader<Grid>::read(factory, mshfile, boundary, *physical, true, false);

        // Create the grid
        grid = std::shared_ptr<Grid>(factory.createGrid());

        // Loadbalance the grid distributing physical entity information
        using PhysicalMap = std::map<typename Grid::GlobalIdSet::IdType, int>;
        const auto& idset = grid->globalIdSet();
        const auto& indexset = grid->leafGridView().indexSet();
        PhysicalMap id_to_data;
        for(auto e : elements(grid->leafGridView()))
          id_to_data[idset.id(e)] = (*physical)[indexset.index(e)];

        impl::PhysicalInfoDataHandle datahandle(idset, id_to_data);
        grid->loadBalance(datahandle);

        // Undo the mapping nature of the physical data
        const auto& newindexset = grid->leafGridView().indexSet();
        physical->resize(grid->size(0));
        for(auto e : elements(grid->leafGridView()))
          (*physical)[newindexset.index(e)] = id_to_data[idset.id(e)];
      }

      return grid;
    }

    Parameter createParameters()
    {
      // If createParameters is called before createGrid, we make sure
      // to call createGrid from here!
      if (!physical)
	createGrid();

      return std::make_tuple(physical);
    }

    private:
    Dune::ParameterTree config;
    std::shared_ptr<Grid> grid;
    std::shared_ptr<std::vector<int>> physical;
  };

} // namespace Dune::BlockLab

#endif
