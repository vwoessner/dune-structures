#ifndef DUNE_STRUCTURES_PRESTRESS_HH
#define DUNE_STRUCTURES_PRESTRESS_HH

#include<dune/common/fmatrix.hh>
#include<dune/common/parametertree.hh>

#include<memory>
#include<numeric>
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
class StraightFibrePrestress
  : public MaterialPrestressBase<GV, T>
{
  public:
  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;

  StraightFibrePrestress(const Dune::ParameterTree& param)
    : scale(param.get<T>("scale"))
    , dir(param.get<Dune::FieldVector<T, 3>>("direction"))
  {
    T norm = dir.two_norm();
    for (int i=0; i<3; ++i)
      dir[i] = dir[i] / norm;
  }

  virtual ~StraightFibrePrestress() {}

  virtual void evaluate(const Entity&, const Coord&, Dune::FieldMatrix<T, 3, 3>& m) const override
  {
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
        m[i][j] = scale * dir[i] * dir[j];
  }

  private:
  T scale;
  Dune::FieldVector<T, 3> dir;
};


template<typename GV, typename T>
class CurvedFibrePrestress
  : public MaterialPrestressBase<GV, T>
{
  public:
  using Entity = typename GV::template Codim<0>::Entity;
  using Coord = typename GV::template Codim<0>::Geometry::LocalCoordinate;

  CurvedFibrePrestress(const Dune::ParameterTree& param, const Dune::ParameterTree& rootconfig)
    : scale(param.get<T>("scale", T(1.0)))
    , sampling(param.get<int>("sampling", 50))
    , subsampling(param.get<int>("subsampling", 100))
  {
    if (!param.hasKey("fibre"))
      DUNE_THROW(MaterialError, "CurvedFibrePrestress expects a fibre config!");

    auto& fibreconfig = rootconfig.sub(param.get<std::string>("fibre"));
    if (fibreconfig.get<std::string>("shape") != "overnucleus")
      DUNE_THROW(MaterialError, "CurvedFibrePrestress only supports overnucleus fibres!");

    // Extracting control points from config. This makes quite some assumptions of course.
    auto scale = rootconfig.get<T>("grid.scaling", T(1.0));
    auto start = fibreconfig.get<Dune::FieldVector<T, 3>>("start") * scale;
    auto middle = fibreconfig.get<Dune::FieldVector<T, 3>>("middle") * scale;
    auto end = fibreconfig.get<Dune::FieldVector<T, 3>>("end") * scale;
    auto slope = fibreconfig.get<T>("slope");

    // Resize the control point data structure
    control_points.resize(2);
    control_points[0].resize(3);
    control_points[1].resize(3);

    // Add control points
    control_points[0][0] = start;
    control_points[0][1] = Dune::FieldVector<T, 3>{{(1.0 - slope) * start[0] + slope * middle[0],
                                                    (1.0 - slope) * start[1] + slope * middle[1],
                                                    middle[2]}};
    control_points[0][2] = middle;

    control_points[1][0] = middle;
    control_points[1][1] = Dune::FieldVector<T, 3>{{(1.0 - slope) * middle[0] + slope * end[0],
                                                    (1.0 - slope) * middle[1] + slope * end[1],
                                                    middle[2]}};
    control_points[1][2] = end;
  }

  virtual ~CurvedFibrePrestress() {}

  virtual void evaluate(const Entity& e, const Coord& x, Dune::FieldMatrix<T, 3, 3>& m) const override
  {
    auto xg = e.geometry().global(x);
    auto [cpindex, pos] = closest_point(xg);
    auto dir = bezier_curve_tangent(cpindex, pos);

    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
        m[i][j] = scale * dir[i] * dir[j];
  }

  // This function can be used for debugging purposes. It calculates the minimum
  // distance point on the curve and then returns the distance to that point.
  // Visualized, this can give insight into what we are doing.
  T distance_to_minimum(const Entity& e, const Coord& x) const
  {
    auto xg = e.geometry().global(x);
    auto [cpindex, pos] = closest_point(xg);
    return distance(xg, cpindex, pos);
  }

  private:
  T bernstein_polynomial(int n, int i, T x) const
  {
    // Necessary for the derivative, not caring about performance here
    if ((i<0) || (i>n))
      return 0.0;

    // Initialize an empty product
    T prod(1.0);

    // Multiply with n choose i
    for (int j=1; j<=i; ++j)
      prod *= T(n - (i - j)) / T(j);

    // Multiply with x^i * (1-x)^(n-i)
    using std::pow;
    return prod * std::pow(x, i) * std::pow(1.0 - x, n-i);
  }

  T bernstein_polynomial_diff(int n, int i, T x) const
  {
    return n * (bernstein_polynomial(n-1, i-1, x) - bernstein_polynomial(n-1, i, x));
  }

  Dune::FieldVector<T, 3> bezier_curve(int cpindex, T x) const
  {
    Dune::FieldVector<T, 3> ret;
    int n = control_points[cpindex].size() - 1;
    for(int i=0; i <= n; ++i)
      ret.axpy(bernstein_polynomial(n, i, x), control_points[cpindex][i]);
    return ret;
  }

  Dune::FieldVector<T, 3> bezier_curve_tangent(int cpindex, T x) const
  {
    Dune::FieldVector<T, 3> ret;
    int n = control_points[cpindex].size() - 1;
    for(int i=0; i <= n; ++i)
      ret.axpy(bernstein_polynomial_diff(n, i, x), control_points[cpindex][i]);

    auto norm = ret.two_norm();
    ret /= norm;
    return ret;
  }


  T distance(const Dune::FieldVector<T, 3>& point, int cpindex, T x) const
  {
    auto curve = bezier_curve(cpindex, x);
    curve -= point;
    return curve.two_norm();
  }

  T minimize_across_curve(const Dune::FieldVector<T, 3>& point, int cpindex) const
  {
    T h = 1.0 / sampling;
    std::vector<T> sample_points(sampling);
    for (int i=0; i<=sampling; ++i)
      sample_points[i] = i * h;

    using std::numeric_limits;
    auto cmp = [this, point, cpindex](auto a, auto b) {
      return this->distance(point, cpindex, a) < this->distance(point, cpindex, b);
    };

    T close = *std::min_element(sample_points.begin(),
                                sample_points.end(),
                                cmp);

    auto hsub = 2.0 * (h / subsampling);
    sample_points.resize(subsampling);
    for (int i=0; i<=subsampling; ++i)
      sample_points[i] = close - h + i * hsub;

    return *std::min_element(sample_points.begin(),
                             sample_points.end(),
                             cmp);
  }

  std::pair<int, T> closest_point(const Dune::FieldVector<T, 3>& point) const
  {
    int curve = 0;
    T minpos = minimize_across_curve(point, 0);
    T mindist = distance(point, 0, minpos);

    for(int i=1; i<control_points.size(); ++i)
    {
      auto pos = minimize_across_curve(point, i);
      auto dist = distance(point, i, pos);
      if (dist < mindist)
      {
        curve = i;
        minpos = pos;
        mindist = dist;
      }
    }

    return std::make_pair(curve, minpos);
  }

  private:
  std::vector<std::vector<Dune::FieldVector<T, 3>>> control_points;
  T scale;
  int sampling, subsampling;
};


template<typename GV, typename T>
std::shared_ptr<MaterialPrestressBase<GV, T>> construct_prestress(const Dune::ParameterTree& params, const Dune::ParameterTree& rootparams)
{
  auto type = params.get<std::string>("type", "none");

  if (type == "isotropic")
    return std::make_shared<IsotropicPrestress<GV, T>>(params);

  if (type == "straight_fibre")
    return std::make_shared<StraightFibrePrestress<GV, T>>(params);

  if (type == "curved_fibre")
    return std::make_shared<CurvedFibrePrestress<GV, T>>(params, rootparams);

  if (type != "none")
    DUNE_THROW(MaterialError, "Your material specification is not supported!");

  return std::make_shared<NoPrestress<GV, T>>();
}

#endif
