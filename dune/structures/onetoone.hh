#ifndef DUNE_STRUCTURES_ONETONE_HH
#define DUNE_STRUCTURES_ONETONE_HH

/** Tools to check that a given displacement field maps the reference
 *  configuration in a one-to-one fashion. Failure to do so results in
 *  unphysical solutions that cause all sorts of numerical trouble. This
 *  is intended as a debugging tool to diagnose such failures.
 */

#include<dune/common/fvector.hh>
#include<dune/pdelab.hh>

#include<array>


bool below_plane(std::array<Dune::FieldVector<double, 3>, 4> points, int index)
{
  auto& p = points[index];
  auto& x = points[(index < 1) ? 1 : 0];
  auto& y = points[(index < 2) ? 2 : 1];
  auto& z = points[(index < 3) ? 3 : 2];

  Dune::FieldVector<double, 3> n;
  Dune::PDELab::crossproduct(y-x, z-x, n);

  double zeval = (n * x - n[0] * p[0] - n[1] * p[1]) / n[2];

  return p[2] < zeval;
}


bool fail_check(std::array<Dune::FieldVector<double, 3>, 4> points1,
                std::array<Dune::FieldVector<double, 3>, 4> points2)
{
  bool failed = false;

  failed = failed || (below_plane(points1, 0) != below_plane(points2, 0));
  failed = failed || (below_plane(points1, 1) != below_plane(points2, 1));
  failed = failed || (below_plane(points1, 2) != below_plane(points2, 2));
  failed = failed || (below_plane(points1, 3) != below_plane(points2, 3));

  return failed;
}


/** Checks whether the displacement field defined by the given vector grid function
 *  is one-to-one. Restricted to P1 continuous finite elements.
 */
template<typename GF>
bool is_onetoone(const GF& gf, bool verbose=false)
{
  for (const auto& e: elements(gf.getGridView()))
  {
    auto geo = e.geometry();
    auto ref = Dune::referenceElement(geo);

    std::array<Dune::FieldVector<double, 3>, 4> ref_corners;
    std::array<Dune::FieldVector<double, 3>, 4> displ_corners;

    for (int i=0; i<4; ++i)
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
        for (int i=0; i<4; ++i)
          std::cout << "Corner " << i << ": " << ref_corners[i] << " -> " << displ_corners[i] << std::endl;
        std::cout << std::endl;
      }
      return false;
    }
  }

  return true;
}

#endif
