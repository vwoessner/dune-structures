#include "config.h"

#include <exception>
#include <memory>
#include <regex>
#include <string>
#include <tuple>
#include <vector>

#include <yaml-cpp/yaml.h>

#include <dune/common/exceptions.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/blocklab/construction/context.hh>
#include <dune/blocklab/grids.hh>
#include <dune/blocklab/init.hh>
#include <dune/blocklab/vectors.hh>

#include <dune/structures/material.hh>

int
main(int argc, char** argv)
{
  auto init = Dune::BlockLab::initBlockLab(argc, argv);
  Dune::TestSuite test;

  // Valid parameterization
  std::string config_str = "solver:\n"
                           "  grid:\n"
                           "    refinement: 1\n"
                           "    _blockname: grid_0\n"
                           "  blocks:\n"
                           "    material:\n"
                           "      _blockname: material_0\n"
                           "      debug: True\n"
                           "      materials:\n"
                           "        -\n"
                           "          group: 0\n"
                           "          model: linear\n"
                           "          poisson_ratio: 0.3333\n"
                           "          youngs_modulus: 300\n";
  YAML::Node config;
  config = YAML::Load(config_str);

  // Grid and Vector
  using Grid = Dune::UGGrid<2>;
  using GridProvider = Dune::BlockLab::StructuredSimplexGridProvider<Grid>;
  using VectorProvider = Dune::BlockLab::PkFemVectorProvider<GridProvider, 1>;
  auto grid = std::make_shared<GridProvider>(config);
  auto vector = std::make_shared<VectorProvider>(grid);

  // Construction Context
  using Material = std::shared_ptr<
    ElasticMaterialBase<Dune::PDELab::OverlappingEntitySet<
                          typename GridProvider::Grid::Traits::LeafGridView>,
                        double>>;
  using ParameterTuple =
    std::tuple<std::shared_ptr<std::vector<int>>, Material>;

  using Context =
    Dune::BlockLab::ConstructionContext<ParameterTuple, VectorProvider>;

  // Valid parameterization
  try
  {
    Context ctx(init.helper, config, vector);
    ctx.template registerBlock<MaterialInitializationBlock>("material");
    auto solver = ctx.constructSolver(config["solver"]);
    solver->apply();
    test.check(true, "Valid Parameterization");
  }
  catch (std::exception& e)
  {
    test.check(false, "Valid Parameterization")
      << "Valid parameterization could not be loaded: " << e.what();
  }

  // Replace group by "1", invalidating the parameterization
  std::regex group_re("group: (\\d+)");
  std::smatch group_match;
  std::regex_search(config_str, group_match, group_re);
  config_str.replace(group_match[1].first, group_match[1].second, "1", 1);
  config = YAML::Load(config_str);

  // Invalid parameterization
  try
  {
    Context ctx(init.helper, config, vector);
    ctx.template registerBlock<MaterialInitializationBlock>("material");
    auto solver = ctx.constructSolver(config["solver"]);
    solver->apply();
    test.check(false, "Invalid Parameterization Group")
      << "Invalid parameterization was loaded successfully!";
  }
  catch (Dune::IOError& e)
  {
    const std::string error_msg(e.what());
    const auto pos = error_msg.find(
      "Material group 0 is missing in the material configuration!");
    test.check(pos != std::string::npos, "Invalid Parameterization Group")
      << "Invalid parameterization led to wrong error message: " << e.what();
  }
  catch (...)
  {
    test.check(false, "Invalid Parameterization Group")
      << "Invalid parameterization led to wrong error type!";
  }

  // Replace debug with "False"
  std::regex debug_re("debug: (\\w+)");
  std::smatch debug_match;
  std::regex_search(config_str, debug_match, debug_re);
  config_str.replace(debug_match[1].first, debug_match[1].second, "False");
  config = YAML::Load(config_str);

  // Invalid parameterization without debugging
  try
  {
    Context ctx(init.helper, config, vector);
    ctx.template registerBlock<MaterialInitializationBlock>("material");
    auto solver = ctx.constructSolver(config["solver"]);
    solver->apply();
    test.check(true, "Invalid Parameterization Group, Debug off")
      << "Invalid parameterization was not loaded successfully, although debug "
         "was disabled!";
  }
  catch (...)
  {
    test.check(false, "Invalid Parameterization Group, Debug off");
  }

  return test.exit();
}
