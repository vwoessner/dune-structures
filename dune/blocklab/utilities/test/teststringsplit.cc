#include"config.h"

#include<dune/blocklab/utilities/stringsplit.hh>

#include<iostream>
#include<string>
#include<vector>


int check(std::string str, std::vector<std::string> ref)
{
  auto split = Dune::BlockLab::string_split(str);

  if (split.size() != ref.size())
  {
    std::cout << "Split of '" << str << "' did yield " << split.size() << " tokes, expected " << ref.size() << std::endl;
    return 1;
  }

  int count = 0;
  for(std::size_t i=0; i<ref.size(); ++i)
    if (split[i] != ref[i])
      ++count;

  if(count)
  {
    std::cout << "Failed to split '" << str << "':" << std::endl;
    std::cout << "  Got:      ";
    for (auto s : split)
      std::cout << "'" << s << "' ";
    std::cout << std::endl;
    std::cout << "  Expected: ";
    for (auto s : ref)
      std::cout << "'" << s << "' ";
    std::cout << std::endl;
  }

  return count;
}

int main()
{
  int failures = 0;

  failures += check("  a, b,c ", {"a","b","c"});

  return failures;
}
