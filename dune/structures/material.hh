#ifndef DUNE_STRUCTURES_MATERIAL_HH
#define DUNE_STRUCTURES_MATERIAL_HH

#include<dune/common/shared_ptr.hh>
#include<dune/common/parametertree.hh>

#include<map>
#include<memory>
#include<vector>


template<typename GV, typename T>
class ElasticMaterialBase
{
  public:
  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;

  virtual ~ElasticMaterialBase() {}

  virtual T first_lame(const Entity& e, const Coord& x) const = 0;
  virtual T second_lame(const Entity& e, const Coord& x) const = 0;
};


template<typename GV, typename T>
class HomogeneousElasticMaterial : public ElasticMaterialBase<GV, T>
{
  public:
  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;

  // Construct from a parameter tree
  HomogeneousElasticMaterial(const Dune::ParameterTree& params)
	: lame1(params.get<T>("lame1", T(1.0))),
	  lame2(params.get<T>("lame2", T(1.0)))
  {}

  // Construct given expicit parameters
  HomogeneousElasticMaterial(T lame1, T lame2)
    : lame1(lame1),
	  lame2(lame2)
  {}

  virtual T first_lame(const Entity& e, const Coord& x) const override
  {
    return lame1;
  }

  virtual T second_lame(const Entity& e, const Coord& x) const override
  {
    return lame2;
  }

  private:
  T lame1;
  T lame2;
};


template<typename GV, typename T>
class MaterialCollection : public ElasticMaterialBase<GV, T>
{
  public:
  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;

  MaterialCollection(const GV& gv, std::shared_ptr<std::vector<int>> physical_groups)
    : is(gv.indexSet()), physical_entity_mapping(physical_groups)
  {}

  MaterialCollection(const GV& gv, std::vector<int>& physical_groups)
    : is(gv.indexSet()), physical_entity_mapping(Dune::stackobject_to_shared_ptr(physical_groups))
  {}

  void add_material(int material_index, std::shared_ptr<ElasticMaterialBase<GV, T>> material)
  {
    materials.insert(std::make_pair(material_index, material));
  }

  void add_material(int material_index, ElasticMaterialBase<GV, T>& material)
  {
    add_material(material_index, Dune::stackobject_to_shared_ptr(material));
  }

  virtual T first_lame(const Entity& e, const Coord& x) const override
  {
    return get_material(e)->first_lame(e, x);
  }

  virtual T second_lame(const Entity& e, const Coord& x) const override
  {
    return get_material(e)->second_lame(e, x);
  }

  private:
  std::shared_ptr<ElasticMaterialBase<GV, T>> get_material(const Entity& e) const
  {
    return materials.find((*physical_entity_mapping)[is.index(e)])->second;
  }

  const typename GV::IndexSet& is;
  std::shared_ptr<std::vector<int>> physical_entity_mapping;
  std::map<int, std::shared_ptr<ElasticMaterialBase<GV, T>>> materials;
};


#endif
