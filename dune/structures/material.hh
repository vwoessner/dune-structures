#ifndef DUNE_STRUCTURES_MATERIAL_HH
#define DUNE_STRUCTURES_MATERIAL_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/utilities/stringsplit.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/fvector.hh>
#include<dune/common/shared_ptr.hh>
#include<dune/common/parametertree.hh>
#include<dune/structures/prestress.hh>

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
T lookup_with_conversion(const YAML::Node& params, std::string name)
{
  if (params[name])
    return params[name].as<T>();

  T lame1, lame2;
  if (params["youngs_modulus"] && params["poisson_ratio"])
  {
    T young = params["youngs_modulus"].as<T>();
    T pr = params["poisson_ratio"].as<T>();

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
  static constexpr int dim = GV::dimension;

  ElasticMaterialBase(const GV& gv) : gv(gv)
  {}

  virtual ~ElasticMaterialBase() {}

  virtual T parameter(const Entity& e, const Coord& x, int i) const = 0;

  virtual void prestress(const Entity& e, const Coord& x, Dune::FieldMatrix<T, dim, dim>&) const = 0;

  virtual int material_law_index(const Entity& e) const = 0;

  T parameter_unrolled(const Entity& e, int i, T x...) const
  {
    return this->parameter(e, Dune::FieldVector<T, dim>{x}, i);
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
  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;
  static constexpr int dim = GV::dimension;

  virtual ~HomogeneousElasticMaterial() {}

  // Construct from a parameter tree
  HomogeneousElasticMaterial(const GV& gv,
                             const YAML::Node& params,
                             const YAML::Node& rootparams)
    : ElasticMaterialBase<GV, T>(gv)
    , prestr(construct_prestress<GV, T>(params["prestress"], rootparams))
  {
    auto lawstr = params["model"].as<std::string>("linear");
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

  virtual void prestress(const Entity& e, const Coord& x, Dune::FieldMatrix<T, dim, dim>& m) const override
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
  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;
  static constexpr int dim = GV::dimension;

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

  virtual void prestress(const Entity& e, const Coord& x, Dune::FieldMatrix<T, dim, dim>& m) const override
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
    // Find the coarse level entity that can be used to look up the physical entity index
    auto lookup_entity = e;
    while (lookup_entity.hasFather())
      lookup_entity = lookup_entity.father();

    return materials.find((*physical_entity_mapping)[is->index(lookup_entity)])->second;
  }

  const typename GV::IndexSet* is;
  std::shared_ptr<std::vector<int>> physical_entity_mapping;
  std::map<int, std::shared_ptr<ElasticMaterialBase<GV, T>>> materials;
};


template<typename T,typename GV>
std::shared_ptr<ElasticMaterialBase<GV, T>> parse_material(
    const GV& gv,
    std::shared_ptr<std::vector<int>> physical_groups,
    const YAML::Node& params,
    const YAML::Node& rootparams
    )
{
  auto coll = std::make_shared<MaterialCollection<GV, T>>(gv, physical_groups);

  for (auto matconfig : params)
  {
    auto material = std::make_shared<HomogeneousElasticMaterial<GV, T>>(gv, matconfig, rootparams);
    coll->add_material(matconfig["group"].as<int>(), material);
  }

  return coll;
}


template<typename P, typename V>
class MaterialInitializationBlock
  : public Dune::BlockLab::BlockBase<P, V>
{
  public:
  using Traits = Dune::BlockLab::BlockTraits<P, V>;
  using Material = std::shared_ptr<ElasticMaterialBase<typename Traits::EntitySet, double>>;

  template<typename Context>
  MaterialInitializationBlock(Context& ctx, const YAML::Node& config)
    : Dune::BlockLab::BlockBase<P, V>(ctx, config)
    , root_config(ctx.getRootConfig())
    , material_config(config["materials"])
  {}

  virtual ~MaterialInitializationBlock() = default;

  virtual void setup() override
  {
    // Make sure that there is a physical vector - even if the grid provider does not have it
    auto size = this->solver->template getVector<0>()->gridFunctionSpace().gridView().size(0);
    this->solver->introduce_parameter("physical", std::make_shared<std::vector<int>>(size, 0));

    // And then intialize the material
    material = parse_material<double>(
      this->solver->template getVector<0>()->gridFunctionSpace().entitySet(),
      this->solver->template param<std::shared_ptr<std::vector<int>>>("physical"),
      material_config,
      root_config
    );
    this->solver->introduce_parameter("material", typename Traits::Parameter(material));
  }

  static std::vector<std::string> blockData()
  {
    auto data = Dune::BlockLab::BlockBase<P, V>::blockData();
    data.push_back(
      "title: Elasticity Material Initialization                   \n"
      "category: structures                                        \n"
      "schema:                                                     \n"
      "  materials:                                                \n"
      "    type: list                                              \n"
      "    schema:                                                 \n"
      "      type: dict                                            \n"
      "      schema:                                               \n"
      "        model:                                              \n"
      "          type: string                                      \n"
      "          default: linear                                   \n"
      "          allowed:                                          \n"
      "            - linear                                        \n"
      "          meta:                                             \n"
      "            title: Material Law                             \n"
      "        youngs_modulus:                                     \n"
      "          type: float                                       \n"
      "          default: 1e5                                      \n"
      "          dependencies:                                     \n"
      "            model: linear                                   \n"
      "          meta:                                             \n"
      "            title: Youngs's modulus                         \n"
      "        poisson_ratio:                                      \n"
      "          type: float                                       \n"
      "          default: 0.4                                      \n"
      "          dependencies:                                     \n"
      "            model: linear                                   \n"
      "          meta:                                             \n"
      "            title: Poisson ratio                            \n"
      "        group:                                              \n"
      "          type: integer                                     \n"
      "          default: 0                                        \n"
      "          meta:                                             \n"
      "            title: Physical group                           \n"
      "    meta:                                                   \n"
      "      title: Materials                                      \n"
    );
    return data;
  }
  

  protected:
  YAML::Node root_config;
  YAML::Node material_config;
  Material material;
};


#endif
