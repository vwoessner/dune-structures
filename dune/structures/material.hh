#ifndef DUNE_STRUCTURES_MATERIAL_HH
#define DUNE_STRUCTURES_MATERIAL_HH

#include<dune/common/fmatrix.hh>
#include<dune/common/fvector.hh>
#include<dune/common/shared_ptr.hh>
#include<dune/common/parametertree.hh>
#include<dune/structures/prestress.hh>
#include<dune/structures/utilities.hh>

#include<map>
#include<memory>
#include<string>
#include<vector>


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

/* Same applies here. Actually it is not used on the Python side at all right now */
std::map<std::string, std::map<int, std::string>> param_to_index = {
    {"linear", {{0, "first_lame"},
                {1, "second_lame"}
                }
    },
    {"stvenantkirchhoff", {{0, "first_lame"},
                           {1, "second_lame"}
                           }
    },
    {"neohookean", {{0, "first_lame"},
                    {1, "second_lame"}
                    }
    },
};


template<typename T>
T lookup_with_conversion(const Dune::ParameterTree& params, std::string name)
{
  if (params.hasKey(name))
    return params.get<T>(name);

  T lame1, lame2;
  if (params.hasKey("youngs_modulus") && params.hasKey("poisson_ratio"))
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

  if (name == "first_lame")
    return lame1;
  if (name == "second_lame")
    return lame2;

  DUNE_THROW(MaterialError, "Your material specification is not supported!");
  return T(); // Silences warning
}


template<typename GV, typename T>
class ElasticMaterialBase
{
  public:
  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;

  ElasticMaterialBase(const GV& gv) : gv(gv)
  {}

  virtual ~ElasticMaterialBase() {}

  virtual T parameter(const Entity& e, const Coord& x, int i) const = 0;

  virtual void prestress(const Entity& e, const Coord& x, Dune::FieldMatrix<T, 3, 3>&) const = 0;

  virtual int material_law_index(const Entity& e) const = 0;

  virtual T parameter(const Entity& e, int i, T x0, T x1, T x2) const
  {
    return this->parameter(e, Dune::FieldVector<T, 3>{x0, x1, x2}, i);
  }

  GV gridView() const
  {
    return gv;
  }

  virtual std::shared_ptr<std::vector<int>> getPhysicalGroups() const
  {
    return std::make_shared<std::vector<int>>();
  }

  private:
  GV gv;
};


template<typename GV, typename T>
class HomogeneousElasticMaterial : public ElasticMaterialBase<GV, T>
{
  public:
  // I need this weirdo here for the moment.
  // Godbolt experiment: https://godbolt.org/z/xBB7zC
  using ElasticMaterialBase<GV, T>::parameter;

  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;

  virtual ~HomogeneousElasticMaterial() {}

  // Construct from a parameter tree
  HomogeneousElasticMaterial(const GV& gv,
                             const Dune::ParameterTree& params,
                             const Dune::ParameterTree& rootparams)
    : ElasticMaterialBase<GV, T>(gv)
    , prestr(construct_prestress<GV, T>(params.sub("prestress"), rootparams))
  {
    auto lawstr = params.get<std::string>("model", "linear");
    law = law_to_index[lawstr];

    auto& paramnamemap = param_to_index[lawstr];
    parameters.resize(paramnamemap.size());

    for (auto [index, name] : paramnamemap)
      parameters[index] = lookup_with_conversion<T>(params, name);
  }

  virtual T parameter(const Entity& e, const Coord& x, int i) const override
  {
    return parameters[i];
  }

  virtual void prestress(const Entity& e, const Coord& x, Dune::FieldMatrix<T, 3, 3>& m) const override
  {
    prestr->evaluate(e, x, m);
  }

  virtual int material_law_index(const Entity& e) const override
  {
    return law;
  }

  private:
  std::vector<T> parameters;
  std::shared_ptr<MaterialPrestressBase<GV, T>> prestr;
  int law;
};


template<typename GV, typename T>
class MaterialCollection : public ElasticMaterialBase<GV, T>
{
  public:
  // I need this weirdo here for the moment.
  // Godbolt experiment: https://godbolt.org/z/xBB7zC
  using ElasticMaterialBase<GV, T>::parameter;

  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;

  virtual ~MaterialCollection() {}

  MaterialCollection(const GV& gv, std::shared_ptr<std::vector<int>> physical_groups)
    : ElasticMaterialBase<GV, T>(gv),
      is(&gv.indexSet()),
      physical_entity_mapping(physical_groups)
  {}

  MaterialCollection(const GV& gv, std::vector<int>& physical_groups)
    : ElasticMaterialBase<GV, T>(gv),
      is(&gv.indexSet()),
      physical_entity_mapping(Dune::stackobject_to_shared_ptr(physical_groups))
  {}

  void add_material(int material_index, std::shared_ptr<ElasticMaterialBase<GV, T>> material)
  {
    materials.insert(std::make_pair(material_index, material));
  }

  void add_material(int material_index, ElasticMaterialBase<GV, T>& material)
  {
    add_material(material_index, Dune::stackobject_to_shared_ptr(material));
  }

  virtual T parameter(const Entity& e, const Coord& x, int i) const override
  {
    return get_material(e)->parameter(e, x ,i);
  }

  virtual void prestress(const Entity& e, const Coord& x, Dune::FieldMatrix<T, 3, 3>& m) const override
  {
    get_material(e)->prestress(e, x, m);
  }

  virtual int material_law_index(const Entity& e) const override
  {
    return get_material(e)->material_law_index(e);
  }

  virtual std::shared_ptr<std::vector<int>> getPhysicalGroups() const
  {
    return physical_entity_mapping;
  }

  private:
  std::shared_ptr<ElasticMaterialBase<GV, T>> get_material(const Entity& e) const
  {
    return materials.find((*physical_entity_mapping)[is->index(e)])->second;
  }

  const typename GV::IndexSet* is;
  std::shared_ptr<std::vector<int>> physical_entity_mapping;
  std::map<int, std::shared_ptr<ElasticMaterialBase<GV, T>>> materials;
};


template<typename T,typename GV>
std::shared_ptr<MaterialCollection<GV, T>> parse_material(
    const GV& gv,
    std::shared_ptr<std::vector<int>> physical_groups,
    const Dune::ParameterTree& params,
    const Dune::ParameterTree& rootparams
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
    auto material = std::make_shared<HomogeneousElasticMaterial<GV, T>>(gv, groupconf, rootparams);
    coll->add_material(groupconf.get<int>("group", 0), material);
  }

  return coll;
}

#endif
