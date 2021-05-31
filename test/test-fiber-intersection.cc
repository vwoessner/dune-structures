#include "config.h"

#include <exception>
#include <iostream>
#include <memory>

#include <yaml-cpp/yaml.h>

#include <dune/blocklab/grids.hh>
#include <dune/blocklab/init.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/grid/uggrid.hh>

#include <dune/structures/parametrizedcurves.hh>

int
main(int argc, char* argv[])
{
  try
  {
    auto init = Dune::BlockLab::initBlockLab(argc, argv);
    Dune::TestSuite test;
    YAML::Node config;

    // Create the grid
    constexpr int dim = 2;
    using Grid = Dune::UGGrid<dim>;
    using GridProvider = Dune::BlockLab::StructuredSimplexGridProvider<Grid>;
    config["grid"]["N"] = std::vector<int>({ 1, 1 });
    GridProvider gp(config["grid"]);
    const auto grid = gp.createGrid();

    using namespace Dune::FloatCmp;

    // Fiber starting on boundary, ending on intersection
    Dune::TestSuite test_1("Fiber starts at boundary, ends at intersection");
    {
      config["fiber"]["start"] = std::vector<double>({ 0.0, 0.5 });
      config["fiber"]["end"] = std::vector<double>({ 0.5, 0.5 });
      const auto fiber = std::make_shared<StraightFibre<dim>>(config["fiber"]);

      // Run the intersection algorithm
      const auto eps = 1e-12;
      const FibreGridIntersection fiber_data =
        compute_grid_fibre_intersection(grid->leafGridView(), fiber, eps);
      const auto element_is = fiber_data.element_fibre_intersections;
      const auto facet_is = fiber_data.facet_fibre_intersections;

      // Test the result
      test_1.check(element_is.size() == 1, "Number of element intersections");
      test_1.check(facet_is.size() == 1, "Number of facet intersections");

      // Element intersections
      test_1.require(element_is.find(1) != element_is.end(), "Element index");
      const auto [el_start, el_end] = element_is.at(1);
      test_1.check(eq(el_start, 0.0, eps), "Element intersection start");
      test_1.check(eq(el_end, 1.0, eps), "Element intersection end");

      // Facet intersections
      const auto idx = std::make_pair(1, 0);
      test_1.require(facet_is.find(idx) != facet_is.end(),
                     "Facet element indices");
      test_1.check(eq(facet_is.at(idx), 1.0, eps), "Facet intersection");
    }
    test.subTest(test_1);

    // Fiber starting on intersection, ending outside boundary
    Dune::TestSuite test_2(
      "Fiber starts at intersection, ends outside boundary");
    {
      config["fiber"]["start"] = std::vector<double>({ 0.5, 0.5 });
      config["fiber"]["end"] = std::vector<double>({ 1.5, 0.5 });
      const auto fiber = std::make_shared<StraightFibre<dim>>(config["fiber"]);

      // Run the intersection algorithm
      const auto eps = 1e-10;
      const FibreGridIntersection fiber_data =
        compute_grid_fibre_intersection(grid->leafGridView(), fiber, eps);
      const auto element_is = fiber_data.element_fibre_intersections;
      const auto facet_is = fiber_data.facet_fibre_intersections;

      // Test the result
      test_2.check(element_is.size() == 1, "Number of element intersections");
      test_2.check(facet_is.size() == 1, "Number of facet intersections");

      // Element intersections
      test_2.require(element_is.find(0) != element_is.end(), "Element index");
      const auto [el_start, el_end] = element_is.at(0);
      test_2.check(eq(el_start, 0.0, eps), "Element intersection start");
      test_2.check(eq(el_end, 0.5, eps), "Element intersection end");

      // Facet intersections
      const auto idx = std::make_pair(1, 0);
      test_2.require(facet_is.find(idx) != facet_is.end(),
                     "Facet element indices");
      test_2.check(eq(facet_is.at(idx), 0.0, eps), "Facet intersection");
    }
    test.subTest(test_2);

    return test.exit();
  }
  catch (std::exception& e)
  {
    std::cerr << "Error during testing: " << e.what() << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << "Unknown error during testing!" << std::endl;
    return 1;
  }
}
