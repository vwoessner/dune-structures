#ifndef DUNE_BLOCKLAB_CONSTRUCTION_DATAPARSER_HH
#define DUNE_BLOCKLAB_CONSTRUCTION_DATAPARSER_HH

/** Some utilities to parse parameters from strings
 *  into data types that are supported by the parameter
 *  system. These can be used in block construction
 *  from ini files.
 */
#include<dune/common/exceptions.hh>
#include<dune/common/parametertree.hh>
#include<string>

namespace Dune::BlockLab {

  template<typename T>
  T parse_data(const std::string& s)
  {
    return Dune::ParameterTree::Parser<T>::parse(s);
  }

  template<typename P>
  P parse_parameter(const Dune::ParameterTree& param)
  {
    auto type = param.get<std::string>("datatype", "double");
    auto value = param.get<std::string>("value");

    if (type == "double")
      return P(parse_data<double>(value));
    else if (type == "int")
      return P(parse_data<int>(value));
    else if (type == "string")
      return P(value);
    else
      DUNE_THROW(Dune::Exception, "Cannot parse parameter");
  }

} // namespace Dune::BlockLab

#endif
