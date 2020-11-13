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
    ret["filename"] = modconfig["filename"];
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
      "  cytoplasm:                                         \n"
      "    type: dict                                       \n"
      "    schema:                                          \n"
      "      shape:                                         \n"
      "        type: string                                 \n"
      "        allowed:                                     \n"
      "          - box                                      \n"
      "          - sphere                                   \n"
      "          - round                                    \n"
      "          - spread                                   \n"
      "          - ellipsoid                                \n"
      "        default: box                                 \n"
      "        meta:                                        \n"
      "          title: Shape                               \n"
      "      radius:                                        \n"
      "        type: float                                  \n"
      "        default: 1.0                                 \n"
      "        dependencies:                                \n"
      "          shape:                                     \n"
      "            - sphere                                 \n"
      "            - round                                  \n"
      "            - spread                                 \n"
      "        meta:                                        \n"
      "          title: Radius                              \n"
      "      center:                                        \n"
      "        type: list                                   \n"
      "        minlength: " + std::to_string(dim) + "       \n"
      "        maxlength: " + std::to_string(dim) + "       \n"
      "        dependencies:                                \n"
      "          shape:                                     \n"
      "            - sphere                                 \n"
      "            - round                                  \n"
      "            - ellipsoid                              \n"
      "        schema:                                      \n"
      "          type: float                                \n"
      "          default: 0.0                               \n"
      "        meta:                                        \n"
      "          title: Center                              \n"
      "      cutoff:                                        \n"
      "        type: float                                  \n"
      "        default: 0.6                                 \n"
      "        dependencies:                                \n"
      "          shape:                                     \n"
      "            - round                                  \n"
      "            - ellipsoid                              \n"
      "        meta:                                        \n"
      "          title: Sphere Cutoff                       \n"
      "      height:                                        \n"
      "        type: float                                  \n"
      "        default: 1.0                                 \n"
      "        dependencies:                                \n"
      "          shape: spread                              \n"
      "        meta:                                        \n"
      "          title: Height                              \n"
      "      slope:                                         \n"
      "        type: float                                  \n"
      "        default: 0.5                                 \n"
      "        dependencies:                                \n"
      "          shape: spread                              \n"
      "        meta:                                        \n"
      "          title: Slope Parameter                     \n"
      "      radii:                                         \n"
      "        type: list                                   \n"
      "        minlength: " + std::to_string(dim) + "       \n"
      "        maxlength: " + std::to_string(dim) + "       \n"
      "        dependencies:                                \n"
      "          shape: ellipsoid                           \n"
      "        schema:                                      \n"
      "          type: float                                \n"
      "          default: 1.0                               \n"
      "        meta:                                        \n"
      "          title: Ellipsoid Radii                     \n"
      "      lowerleft:                                     \n"
      "        type: list                                   \n"
      "        minlength: " + std::to_string(dim) + "       \n"
      "        maxlength: " + std::to_string(dim) + "       \n"
      "        dependencies:                                \n"
      "          shape: box                                 \n"
      "        schema:                                      \n"
      "          type: float                                \n"
      "          default: 0.0                               \n"
      "        meta:                                        \n"
      "          title: Lowerleft corner                    \n"
      "      upperright:                                    \n"
      "        type: list                                   \n"
      "        minlength: " + std::to_string(dim) + "       \n"
      "        maxlength: " + std::to_string(dim) + "       \n"
      "        dependencies:                                \n"
      "          shape: box                                 \n"
      "        schema:                                      \n"
      "          type: float                                \n"
      "          default: 1.0                               \n"
      "        meta:                                        \n"
      "          title: Upperright corner                   \n"
      "      N:                                             \n"
      "        type: list                                   \n"
      "        minlength: " + std::to_string(dim) + "       \n"
      "        maxlength: " + std::to_string(dim) + "       \n"
      "        dependencies:                                \n"
      "          shape: box                                 \n"
      "        schema:                                      \n"
      "          type: integer                              \n"
      "          default: 10                                \n"
      "          min: 1                                     \n"
      "        meta:                                        \n"
      "          title: Number of grid cells                \n"
      "      group:                                         \n"
      "        type: integer                                \n"
      "        default: 0                                   \n"
      "        min: 0                                       \n"
      "        meta:                                        \n"
      "          title: Physical Group                      \n"
      "      meshwidth:                                     \n"
      "        type: float                                  \n"
      "        default: 0.1                                 \n"
      "        meta:                                        \n"
      "          title: Characteristic Mesh Width           \n"
      "    meta:                                            \n"
      "      title: Cytoplasm                               \n"
      "  scaling:                                           \n"
      "    type: float                                      \n"
      "    default: 1.0                                     \n"
      "    meta:                                            \n"
      "      title: Scale Factor                            \n"
      "  fibres:                                            \n"
      "    type: list                                       \n"
      "    schema:                                          \n"
      "      type: dict                                     \n"
      "      schema:                                        \n"
      "        shape:                                       \n"
      "          type: string                               \n"
      "          allowed:                                   \n"
      "            - cylinder                               \n"
      "            - overnucleus                            \n"
      "          default: cylinder                          \n"
      "          meta:                                      \n"
      "            title: Fiber shape                       \n"
      "        start:                                       \n"
      "          type: list                                 \n"
      "          minlength: " + std::to_string(dim) + "     \n"
      "          maxlength: " + std::to_string(dim) + "     \n"
      "          schema:                                    \n"
      "            type: float                              \n"
      "            default: 0.0                             \n"
      "          meta:                                      \n"
      "            title: Starting Point                    \n"
      "        end:                                         \n"
      "          type: list                                 \n"
      "          minlength: " + std::to_string(dim) + "     \n"
      "          maxlength: " + std::to_string(dim) + "     \n"
      "          schema:                                    \n"
      "            type: float                              \n"
      "            default: 1.0                             \n"
      "          meta:                                      \n"
      "            title: End Point                         \n"
      "        radius:                                      \n"
      "          type: float                                \n"
      "          default: 0.5                               \n"
      "          min: 0                                     \n"
      "          meta:                                      \n"
      "            title: Radius                            \n"
      "        middle:                                      \n"
      "          type: list                                 \n"
      "          minlength: " + std::to_string(dim) + "     \n"
      "          maxlength: " + std::to_string(dim) + "     \n"
      "          dependencies:                              \n"
      "            shape: overnucleus                       \n"
      "          schema:                                    \n"
      "            type: float                              \n"
      "            default: 0.5                             \n"
      "          meta:                                      \n"
      "            title: Mid Point                         \n"
      "        slope:                                       \n"
      "          type: float                                \n"
      "          default: 0.5                               \n"
      "          dependencies:                              \n"
      "            shape: overnucleus                       \n"
      "          meta:                                      \n"
      "            title: Slope Parameter                   \n"
      "        meshwidth:                                   \n"
      "          type: float                                \n"
      "          default: 0.025                             \n"
      "          min: 0                                     \n"
      "          meta:                                      \n"
      "            title: Characteristic Mesh Width         \n"
      "        group:                                       \n"
      "          type: integer                              \n"
      "          default: 1                                 \n"
      "          min: 0                                     \n"
      "          meta:                                      \n"
      "            title: Physical Group                    \n"
      "    meta:                                            \n"
      "      title: Fibers                                  \n"
    };
  }
};

#endif
