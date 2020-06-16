#include"config.h"

#include<dune/blocklab/blocks/filelogger.hh>
#include<dune/blocklab/blocks/parameter.hh>
#include<dune/blocklab/blocks/test/providersetup.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertree.hh>
#include<dune/common/test/testsuite.hh>

#include<fstream>


int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;
  Dune::ParameterTree config;

  config["solver.blocks"] = "parameter, filelogger";
  config["parameter.datatype"] = "int";
  config["parameter.value"] = "42";
  config["parameter.name"] = "myparam";
  config["filelogger.filename"] = "fileloggertest.log";
  config["filelogger.parameter"] = "myparam";
  config["filelogger.datatype"] = "int";

  // Nest this in a block so that the fstream destructor is called
  // for sure before we assert the log file content.
  {
    auto ctx = structured_ug2_p1fem(helper, config);
    ctx.template registerBlock<Dune::BlockLab::ParameterBlock>("parameter");
    ctx.template registerBlock<Dune::BlockLab::FileLoggerBlock>("filelogger");

    auto solver = ctx.constructSolver(config.sub("solver"));
    solver->apply();
  }

  std::fstream stream("fileloggertest.log");
  int parse;
  stream >> parse;

  test.check(parse == 42)
     << "Retrieval of int parameter myparam from log file failed";

  return test.exit();
}
