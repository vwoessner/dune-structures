#include"config.h"

#include<dune/blocklab/blocks/control.hh>
#include<dune/blocklab/blocks/test/providersetup.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertree.hh>
#include<dune/common/test/testsuite.hh>

#include<memory>
#include<tuple>


template<typename P, typename V>
class Dummy
  : public Dune::BlockLab::BlockBase<P, V>
{
  public:
  template<typename Context>
  Dummy(Context&, const Dune::ParameterTree&)
  {}

  virtual ~Dummy() = default;

  virtual void setup() override
  {
    this->solver->template introduce_parameter<int>("count", 0);
  }

  virtual void apply() override
  {
    this->solver->update_parameter("count", this->solver->template param<int>("count") + 1);
  }
};


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;
  Dune::ParameterTree config;

  config["solver.blocks"] = "repeat";
  config["repeat.blocks"] = "dummy";
  config["repeat.iterations"] = "42";

  auto ctx = structured_ug2_p1fem(helper, config);
  ctx.template registerBlock<Dune::BlockLab::RepeatBlock>("repeat");
  ctx.template registerBlock<Dummy>("dummy");

  auto solver = ctx.constructSolver(config.sub("solver"));
  solver->apply();

  test.check(solver->template param<int>("count") == 42)
     << "Dummy block was not executed the correct amount of times";

  return test.exit();
}
