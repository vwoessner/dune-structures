#include"config.h"

#include<dune/blocklab/vectors/concept.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/test/testsuite.hh>

class Vector1
{
  public:
  using FiniteElementMap = int;
};


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;

  test.check(Dune::BlockLab::isVectorProvider<Vector1>())
    << "Vector1 failed the concept check but shouldn't";

  return test.exit();
}
