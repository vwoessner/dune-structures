#ifndef DUNE_BLOCKLAB_CONSTRUCTION_DATAPARSER_HH
#define DUNE_BLOCKLAB_CONSTRUCTION_DATAPARSER_HH

/** Some utilities to parse parameters from strings
 *  into data types that are supported by the parameter
 *  system. These can be used in block construction
 *  from ini files.
 */
#include<dune/blocklab/utilities/stringsplit.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/parametertree.hh>

#include<string>
#include<vector>

namespace Dune::BlockLab {

  template<typename T>
  T parse_data(const std::string& s)
  {
    return Dune::ParameterTree::Parser<T>::parse(s);
  }

  template<typename P>
  P parse_parameter(const std::string& type, const std::string& value)
  {
    if (type == "double")
      return P(parse_data<double>(value));
    else if (type == "int")
      return P(parse_data<int>(value));
    else if (type == "string")
      return P(value);
    else
      DUNE_THROW(Dune::Exception, "Cannot parse parameter");
  }

  template<typename P>
  P parse_parameter(const Dune::ParameterTree& param)
  {
    auto type = param.get<std::string>("datatype", "double");
    auto value = param.get<std::string>("value");

    return parse_parameter<P>(type, value);
  }

  template<typename P>
  std::vector<P> parse_parameter_list(const Dune::ParameterTree& param)
  {
    auto type = param.get<std::string>("datatype", "double");
    auto values = param.get<std::string>("values");

    std::vector<P> ret;
    for (auto value : string_split(values))
      ret.push_back(parse_parameter<P>(type, value));

    return ret;
  }

} // namespace Dune::BlockLab

#endif
