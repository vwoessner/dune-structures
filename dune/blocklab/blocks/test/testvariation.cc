#include"config.h"

#include<dune/blocklab/blocks/variation.hh>
#include<dune/blocklab/blocks/test/providersetup.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertree.hh>
#include<dune/common/test/testsuite.hh>


template<typename P, typename V>
class Dummy
  : public Dune::BlockLab::BlockBase<P, V>
{
  public:
  template<typename Context>
  Dummy(Context&, const Dune::ParameterTree& config)
    : increment(config.get<std::string>("increment"))
  {}

  virtual ~Dummy() = default;

  virtual void setup() override
  {
    this->solver->template introduce_parameter<double>("sum", 0.0);
  }

  virtual void apply() override
  {
    this->solver->update_parameter("sum", this->solver->template param<double>("sum") + this->solver->template param<double>(increment));
  }

  private:
  std::string increment;
};


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;

  {
    Dune::ParameterTree config;
    config["solver.blocks"] = "cv";
    config["cv.blocks"] = "dummy";
    config["cv.iterations"] = "5";
    config["cv.name"] = "foo";
    config["dummy.increment"] = "foo";

    auto ctx = structured_ug2_p1fem(helper, config);
    ctx.template registerBlock<Dummy>("dummy");
    ctx.template registerBlock<Dune::BlockLab::ContinuousVariationBlock>("cv");

    auto solver = ctx.constructSolver(config.sub("solver"));
    solver->apply();

    test.check(std::abs(solver->template param<double>("sum") - 3.0) < 1e-8)
      << "Continuous variation parameter did not sum up as expected!";
  }

  {
    Dune::ParameterTree config;
    config["solver.blocks"] = "dv";
    config["dv.blocks"] = "dummy";
    config["dv.values"] = "0.2, 0.4, 0.6, 0.8, 1.0";
    config["dv.name"] = "foo";
    config["dv.datatype"] = "double";
    config["dummy.increment"] = "foo";

    auto ctx = structured_ug2_p1fem(helper, config);
    ctx.template registerBlock<Dummy>("dummy");
    ctx.template registerBlock<Dune::BlockLab::DiscreteVariationBlock>("dv");

    auto solver = ctx.constructSolver(config.sub("solver"));
    solver->apply();

    test.check(std::abs(solver->template param<double>("sum") - 3.0) < 1e-8)
      << "Discrete variation parameter did not sum up as expected!";
  }

  return test.exit();
}
