#include"config.h"

#include<dune/blocklab/utilities/stringsplit.hh>
#include<dune/common/test/testsuite.hh>

#include<iostream>
#include<string>
#include<vector>


bool check(std::string str, std::vector<std::string> ref)
{
  auto split = Dune::BlockLab::string_split(str);

  if (split.size() != ref.size())
    return false;

  for(std::size_t i=0; i<ref.size(); ++i)
    if (split[i] != ref[i])
      return false;

  return true;
}

int main()
{
  Dune::TestSuite test;

  test.check(check("  a, b,c ", {"a","b","c"})) << "Failed to split string '  a, b,c '";

  return test.exit();
}
