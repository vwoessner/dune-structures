#ifndef DUNE_STRUCTURES_GRIDPROVIDER_HH
#define DUNE_STRUCTURES_GRIDPROVIDER_HH

#include<dune/blocklab/grids/gmsh.hh>
#include<dune/common/parametertree.hh>

#include<cstdlib>
#include<string>


namespace impl {

  template<int dim>
  YAML::Node write_gmsh_config(const YAML::Node& config)
  {
    // Modify the config to also contain a filename etc.
    YAML::Node modconfig(config);
    modconfig["dimension"] = dim;
    modconfig["filename"] = "structures.msh";

    // Write a YML file that only represents the grid
    std::ofstream filestream("gridconfig.yml");
    filestream << config;
    filestream.close();

    // Trigger python mesh generator
    using namespace std::literals;
    auto command = std::string(VIRTUALENV_WRAPPER_DIR) + "/run-in-dune-env generate_cell_mesh gridconfig.yml >meshgen.log"s;
    std::system(command.c_str());

    // And return a config that is understood by dune-blocklab's GMSH grid provider
    YAML::Node ret;
    ret["filename"] = modconfig["filename"].as<std::string>();
    return ret;
  }

} // namespace impl

template<int dim>
class StructuresGridProvider
  : public Dune::BlockLab::GMSHGridProvider<Dune::UGGrid<dim>>
{
  public:
  StructuresGridProvider(const YAML::Node& config)
    : Dune::BlockLab::GMSHGridProvider<Dune::UGGrid<dim>>(impl::write_gmsh_config<dim>(config))
  {}

  static std::vector<std::string> blockData()
  {
    return {
      "title: Structures Cell Grid                          \n"
      "category: grids                                      \n"
      "schema:                                              \n"
      "  cytoplasm_shape:                                   \n"
      "    type: string                                     \n"
      "    allowed:                                         \n"
      "      - box                                          \n"
      "    meta:                                            \n"
      "      title: Cytoplasm Shape                         \n"
      "  lowerleft:                                         \n"
      "    type: list                                       \n"
      "    minlength: " + std::to_string(dim) + "           \n"
      "    maxlength: " + std::to_string(dim) + "           \n"
      "    schema:                                          \n"
      "      type: float                                    \n"
      "      default: 0.0                                   \n"
      "    meta:                                            \n"
      "      title: Lowerleft corner                        \n"
      "  upperright:                                        \n"
      "    type: list                                       \n"
      "    minlength: " + std::to_string(dim) + "           \n"
      "    maxlength: " + std::to_string(dim) + "           \n"
      "    schema:                                          \n"
      "      type: float                                    \n"
      "      default: 1.0                                   \n"
      "    meta:                                            \n"
      "      title: Upperright corner                       \n"
      "  N:                                                 \n"
      "    type: list                                       \n"
      "    minlength: " + std::to_string(dim) + "           \n"
      "    maxlength: " + std::to_string(dim) + "           \n"
      "    schema:                                          \n"
      "      type: integer                                  \n"
      "      default: 10                                    \n"
      "      min: 1                                         \n"
      "    meta:                                            \n"
      "      title: Number of grid cells                    \n"
      "  refinement:                                        \n"
      "    type: integer                                    \n"
      "    default: 0                                       \n"
      "    meta:                                            \n"
      "      title: Global refinement levels                \n"
      "  fibers:                                            \n"
      "    type: list                                       \n"
      "    schema:                                          \n"
      "      type: dict                                     \n"
      "      schema: {}                                     \n"
      "    meta:                                            \n"
      "      title: Fibers                                  \n"
    };
  }
};

#endif
