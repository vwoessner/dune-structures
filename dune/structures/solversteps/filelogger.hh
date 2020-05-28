#ifndef DUNE_STRUCTURES_SOLVERSTEPS_FILELOGGER_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_FILELOGGER_HH

#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/traits.hh>

#include<fstream>


template<typename... V>
class FileLoggerStep
  : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  FileLoggerStep(const Dune::ParameterTree& params)
    : params(params)
  {}

  virtual ~FileLoggerStep() = default;

  virtual void pre() override
  {
    if (params.get<bool>("append", false))
      stream.open(params.get<std::string>("filename"), std::fstream::out | std::fstream::app);
    else
      stream.open(params.get<std::string>("filename"), std::fstream::out);
    if (!stream)
      DUNE_THROW(Dune::Exception, "Cannot open logfile for FileLoggerStep");
  }

  virtual void apply() override
  {
    auto datatype = params.get<std::string>("datatype", "double");
    auto parameter = params.get<std::string>("parameter");

    if (datatype == "double")
      stream << this->solver->template param<double>(parameter) << std::endl;
    else if (datatype == "fieldvector")
      stream << this->solver->template param<Dune::FieldVector<double, Traits::dim>>(parameter) << std::endl;
    else
      DUNE_THROW(Dune::Exception, "Unknown datatype in FileLoggerStep");
  }

  virtual void post() override
  {
    stream.close();
  }

  private:
  Dune::ParameterTree params;
  std::fstream stream;
};

#endif
