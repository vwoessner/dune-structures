#ifndef DUNE_STRUCTURES_PRESTRESS_HH
#define DUNE_STRUCTURES_PRESTRESS_HH

#include<dune/common/fmatrix.hh>
#include<dune/common/parametertree.hh>

#include<memory>
#include<string>


class MaterialError : public Dune::IOError {};


template<typename GV, typename T>
class MaterialPrestressBase
{
  public:
  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;

  virtual ~MaterialPrestressBase() {}

  virtual void evaluate(const Entity&, const Coord&, Dune::FieldMatrix<T, 3, 3>&) const = 0;
};


template<typename GV, typename T>
class NoPrestress
  : public MaterialPrestressBase<GV, T>
{
  public:
  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;

  virtual ~NoPrestress() {}

  virtual void evaluate(const Entity&, const Coord&, Dune::FieldMatrix<T, 3, 3>& m) const override
  {
    m = Dune::FieldMatrix<T, 3, 3>(0.0);
  }
};


template<typename GV, typename T>
class IsotropicPrestress
  : public MaterialPrestressBase<GV, T>
{
  public:
  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;

  IsotropicPrestress(const Dune::ParameterTree& param)
    : scale(param.get<T>("scale"))
  {}

  virtual ~IsotropicPrestress() {}

  virtual void evaluate(const Entity&, const Coord&, Dune::FieldMatrix<T, 3, 3>& m) const override
  {
    Dune::FieldMatrix<T, 3, 3> ret(0.0);
    ret[0][0] = scale;
    ret[1][1] = scale;
    ret[2][2] = scale;
    m = ret;
  }

  private:
  T scale;
};


template<typename GV, typename T>
std::shared_ptr<MaterialPrestressBase<GV, T>> construct_prestress(const Dune::ParameterTree& params)
{
  auto type = params.get<std::string>("type", "none");

  if (type == "isotropic")
    return std::make_shared<IsotropicPrestress<GV, T>>(params);

  if (type != "none")
    DUNE_THROW(MaterialError, "Your material specification is not supported!");

  return std::make_shared<NoPrestress<GV, T>>();
}

#endif
