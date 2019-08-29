#ifndef DUNE_STRUCTURES_GMSH_FACTORY_HH
#define DUNE_STRUCTURES_GMSH_FACTORY_HH

/** A grid factory class for our grids that
 *  * uses the GMSHReader for actually parsing the grid
 *  * distributes physical entity tag information in loadbalancing
 */

#include<dune/common/parametertree.hh>
#include<dune/grid/common/datahandleif.hh>
#include<dune/grid/io/file/gmshreader.hh>

#include<memory>
#include<vector>


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


template<typename Grid>
class PhysicalEntityGmshFactory
{
  public:
  PhysicalEntityGmshFactory(const Dune::MPIHelper& helper, const Dune::ParameterTree& config)
  {
    // Set up data for grid factory
    auto mshfile = config.get<std::string>("gmshfile");
    Dune::GridFactory<Grid> factory;
    std::vector<int> boundary;
    physical = std::make_shared<std::vector<int>>();

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

    PhysicalInfoDataHandle datahandle(idset, id_to_data);
    grid->loadBalance(datahandle);

    // Undo the mapping nature of the physical data
    const auto& newindexset = grid->leafGridView().indexSet();
    physical->resize(grid->size(0));
    for(auto e : elements(grid->leafGridView()))
      (*physical)[newindexset.index(e)] = id_to_data[idset.id(e)];
  }

  std::shared_ptr<Grid> getGrid() const
  {
    return grid;
  }

  std::shared_ptr<std::vector<int>> getPhysical() const
  {
    return physical;
  }
  private:
  std::shared_ptr<Grid> grid;
  std::shared_ptr<std::vector<int>> physical;
};

#endif
