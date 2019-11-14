#ifndef DUNE_STRUCTURES_MATERIAL_HH
#define DUNE_STRUCTURES_MATERIAL_HH

#include<dune/common/fvector.hh>
#include<dune/common/shared_ptr.hh>
#include<dune/common/parametertree.hh>
#include<dune/structures/utilities.hh>

#include<map>
#include<memory>
#include<string>
#include<vector>


class MaterialError : public Dune::IOError {};


/* This mapping has to agree with the mapping in Python.
 * In theory we could code-generate this in order to specify
 * it exactly once, but I do not see much value in that right now.
 */
std::map<std::string, int> law_to_index = {
    {"linear", 0},
    {"stvenantkirchhoff", 1},
    {"neohookean", 2},
    {"mooneyrivlin", 3}
};


template<typename GV, typename T>
class ElasticMaterialBase
{
  public:
  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;

  virtual ~ElasticMaterialBase() {}

  virtual T first_lame(const Entity& e, const Coord& x) const = 0;
  virtual T second_lame(const Entity& e, const Coord& x) const = 0;
  virtual int material_law_index(const Entity& e) const = 0;

  virtual T pretension(const Entity& e, const Coord& x) const
  {
    return T(0.0);
  }

  // These are part of a hack that is soon to go away
  virtual T first_lame(const Entity& e, T x0, T x1, T x2) const
  {
    return this->first_lame(e, Dune::FieldVector<T, 3>{x0, x1, x2});
  }

  virtual T second_lame(const Entity& e, T x0, T x1, T x2) const
  {
    return this->second_lame(e, Dune::FieldVector<T, 3>{x0, x1, x2});
  }

  virtual T pretension(const Entity& e, T x0, T x1, T x2) const
  {
    return this->pretension(e, Dune::FieldVector<T, 3>{x0, x1, x2});
  }
};


template<typename GV, typename T>
class HomogeneousElasticMaterial : public ElasticMaterialBase<GV, T>
{
  public:
  // I need these weirdos here for the moment.
  // Godbolt experiment: https://godbolt.org/z/xBB7zC
  using ElasticMaterialBase<GV, T>::first_lame;
  using ElasticMaterialBase<GV, T>::second_lame;
  using ElasticMaterialBase<GV, T>::pretension;

  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;

  virtual ~HomogeneousElasticMaterial() {}

  // Construct from a parameter tree
  HomogeneousElasticMaterial(const Dune::ParameterTree& params)
  {
    law = law_to_index[params.get<std::string>("model", "linear")];
    pretens = params.get<T>("pretension", T(0.0));
    // The parameters of linear elasticity are ambiguous.
    // We only accept some combinations. They can be taken
    // from here:
    // https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
    if (params.hasKey("lame1") && params.hasKey("lame2"))
    {
      lame1 = params.get<T>("lame1");
      lame2 = params.get<T>("lame2");
    }
    else if (params.hasKey("youngs_modulus") && params.hasKey("poisson_ratio"))
    {
      T young = params.get<T>("youngs_modulus");
      T pr = params.get<T>("poisson_ratio");

      lame1 = (pr * young) / ((1.0 + pr) * (1.0 - 2.0 * pr));
      lame2 = young / (2.0 * (1.0 + pr));
    }
    else
    {
      DUNE_THROW(MaterialError, "Your material specification is not supported!");
    }
  }

  // Construct given explicit Lame parameters
  HomogeneousElasticMaterial(T lame1, T lame2)
    : lame1(lame1),
	  lame2(lame2),
	  law(0)
  {}

  virtual T first_lame(const Entity& e, const Coord& x) const override
  {
    return lame1;
  }

  virtual T second_lame(const Entity& e, const Coord& x) const override
  {
    return lame2;
  }

  virtual T pretension(const Entity& e, const Coord& x) const override
  {
    return pretens;
  }

  virtual int material_law_index(const Entity& e) const override
  {
    return law;
  }

  private:
  T lame1;
  T lame2;
  T pretens;
  int law;
};


template<typename GV, typename T>
class MaterialCollection : public ElasticMaterialBase<GV, T>
{
  public:
  // I need these weirdos here for the moment.
  // Godbolt experiment: https://godbolt.org/z/xBB7zC
  using ElasticMaterialBase<GV, T>::first_lame;
  using ElasticMaterialBase<GV, T>::second_lame;

  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;

  virtual ~MaterialCollection() {}

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

  virtual T pretension(const Entity& e, const Coord& x) const override
  {
    return get_material(e)->pretension(e, x);
  }

  virtual int material_law_index(const Entity& e) const override
  {
    return get_material(e)->material_law_index(e);
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


template<typename T,typename GV>
std::shared_ptr<MaterialCollection<GV, T>> parse_material(
    const GV& gv,
    std::shared_ptr<std::vector<int>> physical_groups,
    const Dune::ParameterTree& params
    )
{
  auto coll = std::make_shared<MaterialCollection<GV, T>>(gv, physical_groups);

  // Get the list of materials and iterate over them
  auto material_groups = params.get<std::string>("materials");
  auto groups = str_split(material_groups);
  for (auto group : groups)
  {
    str_trim(group);
    const auto& groupconf = params.sub(group);
    auto material = std::make_shared<HomogeneousElasticMaterial<GV, T>>(groupconf);
    coll->add_material(groupconf.get<int>("group"), material);
  }

  return coll;
}

#endif
