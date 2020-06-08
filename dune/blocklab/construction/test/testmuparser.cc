#include"config.h"

#include<dune/blocklab/construction/muparser.hh>
#include<dune/blocklab/solver.hh>
#include<dune/common/filledarray.hh>
#include<dune/grid/yaspgrid.hh>

#include<iostream>
#include<memory>
#include<string>


template<typename... Callbacks>
int check(std::string expr, double result, Callbacks&&... callbacks)
{
  constexpr int dim = 3;
  using D = Dune::FieldVector<double, dim>;
  using R = double;

  using Grid = Dune::YaspGrid<dim>;
  D ll(0.0), ur(1.0);
  auto N = Dune::filledArray<dim, unsigned int>(1);
  auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(ll, ur, N);

  auto f = Dune::BlockLab::muparser_callable<R(Grid::LeafGridView::Codim<0>::Entity, D)>(expr, std::forward<Callbacks>(callbacks)...);
  auto eval = f(*grid->leafGridView().template begin<0>(), D(0.5));

  if (std::abs(eval - result) > 1e-8)
  {
    std::cout << "Function '" << expr << "' evaluated to " << eval << "! Expected " << result << std::endl;
    return 1;
  }
  else
    return 0;
}

int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // Construct a solver instance
  auto solver = std::make_shared<Dune::BlockLab::BlockSolver<std::tuple<>, std::tuple<>>>();
  solver->introduce_parameter("foo", 2.0);

  int failures = 0;
  // Check a normal expression that depends on the custom geometric symbols x, y and z
  failures += check("x + y + z", 1.5);

  // Check the function that accesses parameters from the parameter class
  failures += check("param(\"foo\")", 2.0, solver);

  // Check that the evaluation of param(foo) is still correct after changing the value
  // This currently holds, but I am unsure whether I should mark param as volatile in
  // muparser in order to be on the safe side.
  solver->update_parameter("foo", 3.0);
  failures += check("param(\"foo\")", 3.0, solver);

  return 0;
}
