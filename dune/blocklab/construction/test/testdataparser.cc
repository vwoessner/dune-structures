#include"config.h"

#include<dune/blocklab/solver.hh>
#include<dune/blocklab/construction/dataparser.hh>
#include<dune/common/parametertree.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/test/testsuite.hh>


template<typename Checker>
bool check(std::string t, std::string v, Checker&& checker)
{
  using P = Dune::BlockLab::BlockSolver<std::tuple<>, std::tuple<>>::Parameter;

  Dune::ParameterTree config;
  config["datatype"] = t;
  config["value"] = v;

  return checker(Dune::BlockLab::parse_parameter<P>(config));
}

int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;

  test.check(check("double", "42", [](auto x){ return std::abs(std::get<double>(x) - 42.0) < 1e-8; }));
  test.check(check("int", "42", [](auto x){ return std::get<int>(x) == 42; }));

  return test.exit();
}
