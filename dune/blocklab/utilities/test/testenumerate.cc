#include"config.h"

#include<dune/blocklab/utilities/enumerate.hh>

#include<numeric>
#include<vector>


int main()
{
  std::vector<int> v(42);
  std::iota(v.begin(), v.end(), 1);

  int ret = 0;
  for (auto [i, v] : Dune::BlockLab::enumerate(v))
    ret += std::abs(v - (int)i - 1);

  return ret;
}
