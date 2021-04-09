#ifndef DUNE_STRUCTURES_ONETONE_HH
#define DUNE_STRUCTURES_ONETONE_HH

/** Tools to check that a given displacement field maps the reference
 *  configuration in a one-to-one fashion. Failure to do so results in
 *  unphysical solutions that cause all sorts of numerical trouble. This
 *  is intended as a debugging tool to diagnose such failures.
 */

#include <dune/blocklab/blocks/blockbase.hh>
#include <dune/blocklab/blocks/enableif.hh>
#include <dune/common/fvector.hh>
#include <dune/pdelab.hh>

#include <array>

bool
below_plane(std::array<Dune::FieldVector<double, 3>, 4> points, int index)
{
  auto& p = points[index];
  auto& x = points[(index < 1) ? 1 : 0];
  auto& y = points[(index < 2) ? 2 : 1];
  auto& z = points[(index < 3) ? 3 : 2];

  Dune::FieldVector<double, 3> n;
  Dune::PDELab::crossproduct(y - x, z - x, n);

  return n * (p - x) < 0;
}

bool
fail_check(std::array<Dune::FieldVector<double, 3>, 4> points1,
           std::array<Dune::FieldVector<double, 3>, 4> points2)
{
  bool failed = false;

  failed = failed || (below_plane(points1, 0) != below_plane(points2, 0));
  failed = failed || (below_plane(points1, 1) != below_plane(points2, 1));
  failed = failed || (below_plane(points1, 2) != below_plane(points2, 2));
  failed = failed || (below_plane(points1, 3) != below_plane(points2, 3));

  return failed;
}

/** Checks whether the displacement field defined by the given vector grid
 * function is one-to-one. Restricted to P1 continuous finite elements.
 */
template<typename GF>
bool
is_onetoone(const GF& gf, bool verbose = false)
{
  for (const auto& e : elements(gf.getGridView()))
  {
    auto geo = e.geometry();
    auto ref = Dune::referenceElement(geo);

    std::array<Dune::FieldVector<double, 3>, 4> ref_corners;
    std::array<Dune::FieldVector<double, 3>, 4> displ_corners;

    for (int i = 0; i < 4; ++i)
    {
      auto pos = ref.position(i, 3);
      ref_corners[i] = geo.global(pos);
      gf.evaluate(e, pos, displ_corners[i]);
      displ_corners[i] += ref_corners[i];
    }

    if (fail_check(ref_corners, displ_corners))
    {
      if (verbose)
      {
        std::cout << std::endl;
        for (int i = 0; i < 4; ++i)
          std::cout << "Corner " << i << " position: " << ref_corners[i]
                    << " -> " << displ_corners[i] << std::endl;
        std::cout << std::endl;

        for (int i = 0; i < 4; ++i)
          std::cout << "Corner " << i << " check: "
                    << (below_plane(ref_corners, i) ? "below" : "above")
                    << " in ref "
                    << (below_plane(displ_corners, i) ? "below" : "above")
                    << " in displaced " << std::endl;
        std::cout << std::endl;
      }
      return false;
    }
  }

  return true;
}

template<typename P,
         typename V,
         std::size_t i,
         typename enabled = Dune::BlockLab::disabled>
class OneToOneMappingCheckerBlock
  : public Dune::BlockLab::DisabledBlock<P, V, i>
{
public:
  template<typename Context>
  OneToOneMappingCheckerBlock(Context& ctx, const YAML::Node& config)
    : Dune::BlockLab::DisabledBlock<P, V, i>(ctx, config)
  {
  }
};

template<typename P, typename V, std::size_t i>
class OneToOneMappingCheckerBlock<
  P,
  V,
  i,
  Dune::BlockLab::enableBlock<Dune::BlockLab::is3D<P, V, i>()>>
  : public Dune::BlockLab::BlockBase<P, V, i>
{
public:
  using Traits = Dune::BlockLab::BlockTraits<P, V, i>;

  template<typename Context>
  OneToOneMappingCheckerBlock(Context& ctx, const YAML::Node& config)
    : Dune::BlockLab::BlockBase<P, V, i>(ctx, config)
  {
  }

  virtual ~OneToOneMappingCheckerBlock() = default;

  void apply() override
  {
    auto vector = this->solver->template getVector<i>();
    auto& gfs = vector->gridFunctionSpace();
    auto es = gfs.entitySet();

    Dune::PDELab::VectorDiscreteGridFunction vdgf(gfs, *vector);

    std::cout << "Checking one-to-one property of displacement field... ";
    std::cout << (is_onetoone(vdgf) ? "Success!" : "Failure") << std::endl;
  }

  static std::vector<std::string> blockData()
  {
    auto data = Dune::BlockLab::BlockBase<P, V, i>::blockData();
    data.push_back("title: One-to-one displacement mapping checker      \n"
                   "category: structures                                \n");
    return data;
  }
};

#endif
