#include"config.h"

#include<dune/blocklab/utilities/uniquevariant.hh>

int main()
{
  using Variant1 = Dune::BlockLab::unique_variant<int, int>;
  using Variant2 = Dune::BlockLab::unique_variant<double, int, int>;
  using Variant3 = Dune::BlockLab::unique_variant<int, double, int>;

  Variant1 v1(1);
  Variant2 v2(1);
  Variant3 v3(1);

  return 0;
}
