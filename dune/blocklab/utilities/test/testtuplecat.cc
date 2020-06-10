#include"config.h"

#include<dune/blocklab/utilities/tuplecat.hh>

#include<tuple>

int main()
{
  Dune::BlockLab::tuple_cat_t<std::tuple<int, double>, std::tuple<int>> test{1, 1.0, 1};

  return 0;
}
