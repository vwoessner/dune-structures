#ifndef DUNE_STRUCTURES_GRIDPROVIDER_HH
#define DUNE_STRUCTURES_GRIDPROVIDER_HH

#include<dune/blocklab/grids/gmsh.hh>
#include<dune/common/parametertree.hh>

#include<cstdlib>
#include<string>

template<int dimension>
class StructuresGridProvider
  : public Dune::BlockLab::GMSHGridProvider<Dune::UGGrid<dimension>>
{
  public:
  StructuresGridProvider(const YAML::Node& config)
    : Dune::BlockLab::GMSHGridProvider<Dune::UGGrid<dimension>>(config)
  {
    // // Write an inifile that the mesh generator can read
    // Dune::ParameterTree outconfig(config);
    // std::fstream filestream("temporarygrid.ini");
    // config.report(filestream, "grid");
    // filestream.close();

    // // Trigger python mesh generator
    // using namespace std::literals;
    // auto command = std::string(VIRTUALENV_WRAPPER_DIR) + "/run-in-dune-env generate_cell_mesh temporarygrid.ini >meshgen.log"s;
    // std::system(command.c_str());
  }

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
    };
  }
};

#endif
