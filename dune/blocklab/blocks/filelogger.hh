#ifndef DUNE_BLOCKLAB_BLOCKS_FILELOGGER_HH
#define DUNE_BLOCKLAB_BLOCKS_FILELOGGER_HH

/** A block that writes data from the parameter system into
 *  a file. This is a stop-gap measure to collect experiments data.
 *  In the long run, dune-blocklab should incorporate a logging
 *  library.
 */

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/parametertree.hh>

#include<fstream>
#include<string>

namespace Dune::BlockLab {

  template<typename P, typename V>
  class FileLoggerBlock
    : public BlockBase<P, V>
  {
    public:
    using Traits = BlockTraits<P, V>;

    template<typename Context>
    FileLoggerBlock(Context&, const Dune::ParameterTree& config)
      : FileLoggerBlock(config)
    {}

    FileLoggerBlock(const Dune::ParameterTree& config)
      : append(config.get<bool>("append", false))
      , datatype(config.get<std::string>("datatype", "double"))
      , filename(config.get<std::string>("filename"))
      , parameter(config.get<std::string>("parameter"))
    {}

    virtual ~FileLoggerBlock() = default;

    virtual void setup() override
    {
      // Construct the correct flags for the stream
      auto flags = std::fstream::out;
      if (append)
        flags |= std::fstream::app;

      // And open the stream
      stream.open(filename, flags);

      // Throw if opening was not successful
      if (!stream)
        DUNE_THROW(Dune::Exception, "Cannot open logfile for FileLoggerBlock");
    }

    virtual void apply() override
    {
      if (datatype == "double")
        stream << this->solver->template param<double>(parameter) << std::endl;
      else if (datatype == "int")
	stream << this->solver->template param<int>(parameter) << std::endl;
      else
        DUNE_THROW(Dune::Exception, "Unknown datatype '" << datatype << "' in FileLoggerBlock");
    }

    private:
    bool append;
    std::string datatype;
    std::string filename;
    std::string parameter;

    std::fstream stream;
  };

} // namespace Dune::BlockLab

#endif
