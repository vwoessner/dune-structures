#include"config.h"

#include<dune/blocklab/solver.hh>
#include<dune/blocklab/construction/dataparser.hh>
#include<dune/common/parametertree.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/test/testsuite.hh>

#include<numeric>


using P = Dune::BlockLab::BlockSolver<std::tuple<>, std::tuple<>>::Parameter;


template<typename Checker>
bool check(std::string t, std::string v, Checker&& checker)
{
  Dune::ParameterTree config;
  config["datatype"] = t;
  config["value"] = v;

  return checker(Dune::BlockLab::parse_parameter<P>(config));
}

template<typename Checker>
bool check_list(std::string t, std::string v, Checker&& checker)
{
  Dune::ParameterTree config;
  config["datatype"] = t;
  config["values"] = v;

  return checker(Dune::BlockLab::parse_parameter_list<P>(config));
}

int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;

  test.check(
      check("double", "42",
	    [](auto x){
              return std::abs(std::get<double>(x) - 42.0) < 1e-8;
            }));

  test.check(
      check("int", "42",
	    [](auto x){
              return std::get<int>(x) == 42;
            }));

  test.check(check_list("int", "47, 11",
			[](auto vec){
                          return std::get<int>(std::accumulate(vec.begin(), vec.end(), P(0),
						               [](auto a, auto b){
                                                                 return std::get<int>(a) + std::get<int>(b);
                                                               })) == 58;
                        }));

  return test.exit();
}
