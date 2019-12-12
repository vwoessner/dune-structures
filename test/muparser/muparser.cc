#include"config.h"

#include<dune/common/fvector.hh>
#include<dune/structures/muparser.hh>

#include<iostream>


int main()
{
  Dune::FieldVector<double, 3> x{1.0, 2.0, 3.0};
  MuParserCallable c("x*y*z");
  std::cout << c(x) << std::endl;
  return 0;
}
