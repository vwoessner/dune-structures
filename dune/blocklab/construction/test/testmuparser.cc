#include"config.h"

#include<dune/blocklab/construction/muparser.hh>
#include<dune/blocklab/solver.hh>
#include<dune/common/filledarray.hh>
#include<dune/common/test/testsuite.hh>
#include<dune/grid/yaspgrid.hh>

#include<iostream>
#include<memory>
#include<string>


auto make_grid()
{
  constexpr int dim = 3;
  using D = Dune::FieldVector<double, dim>;
  using R = double;

  using Grid = Dune::YaspGrid<dim>;
  D ll(0.0), ur(1.0);
  auto N = Dune::filledArray<dim, unsigned int>(1);
  return Dune::StructuredGridFactory<Grid>::createCubeGrid(ll, ur, N);
}


template<typename... Callbacks>
bool check(std::string expr, double result, Callbacks&&... callbacks)
{
  auto grid = make_grid();
  auto f = Dune::BlockLab::muparser_callable<double(decltype(grid)::element_type::LeafGridView::Codim<0>::Entity, Dune::FieldVector<double, 3>)>(expr, std::forward<Callbacks>(callbacks)...);
  auto eval = f(*grid->leafGridView().template begin<0>(), Dune::FieldVector<double, 3>(0.5));

  return (std::abs(eval - result) < 1e-8);
}


template<typename... Callbacks>
bool check_array(std::string expr, std::array<double, 3> result, Callbacks&&... callbacks)
{
  auto grid = make_grid();
  auto ft = Dune::BlockLab::muparser_callable_array<3, double(decltype(grid)::element_type::LeafGridView::Codim<0>::Entity, Dune::FieldVector<double, 3>)>(expr, std::forward<Callbacks>(callbacks)...);
  std::array<double, 3> error;
  std::transform(ft.begin(), ft.end(), result.begin(), error.begin(),
		 [&grid](auto f, auto res)
		 {
                   return std::abs(f(*grid->leafGridView().template begin<0>(), Dune::FieldVector<double, 3>(0.5)) - res);
		 }
  );

  return (std::accumulate(error.begin(), error.end(), 0.0) < 1e-8);
}

int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;

  // Construct a solver instance
  auto solver = std::make_shared<Dune::BlockLab::BlockSolver<std::tuple<>, std::tuple<>>>(std::tuple<>{}, std::tuple<>{});
  solver->introduce_parameter("foo", 2.0);

  // Check a normal expression that depends on the custom geometric symbols x, y and z
  test.check(check("x + y + z", 1.5)) << "Expression depending on x, y, z failed!";

  // Check the function that accesses parameters from the parameter class
  test.check(check("param(\"foo\")", 2.0, solver)) << "Access to param from solver class failed!";

  // Check that the evaluation of param(foo) is still correct after changing the value
  // This currently holds, but I am unsure whether I should mark param as volatile in
  // muparser in order to be on the safe side.
  solver->update_parameter("foo", 3.0);
  test.check(check("param(\"foo\")", 3.0, solver)) << "Access to param from solver class failed afer update!";;

  // Check that a function array that broadcasts one function to all array entries
  // works correctly
  test.check(check_array("x + y + z", {1.5, 1.5, 1.5}));

  // A function array of three functions
  test.check(check_array("x - 0.5, y, z + 0.5", {0.0, 0.5, 1.0}));

  return 0;
}
