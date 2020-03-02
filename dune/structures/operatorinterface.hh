#ifndef DUNE_STRUCTURES_OPERATORINTERFACE_HH
#define DUNE_STRUCTURES_OPERATORINTERFACE_HH

/** An abstract base class for a local operator.
 *
 * It fixes all the duck-typed arguments in the process.
 */

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/gridfunctionspace/localvector.hh>
#include<dune/pdelab/gridfunctionspace/localfunctionspace.hh>


template<typename GFSU, typename GFSV>
class AbstractLocalOperatorInterface
{
  public:
  using EG = Dune::PDELab::ElementGeometry<typename GFSU::Traits::GridView::template Codim<0>::Entity>;
  using IG = Dune::PDELab::IntersectionGeometry<typename GFSU::Traits::GridView::Intersection>;
  using LFSU = Dune::PDELab::LocalFunctionSpace<GFSU>;
  using LFSV = Dune::PDELab::LocalFunctionSpace<GFSV>;
  using X = Dune::PDELab::LocalVector<double, Dune::PDELab::TrialSpaceTag>;
  using R = Dune::PDELab::WeightedVectorAccumulationView<Dune::PDELab::LocalVector<double, Dune::PDELab::TestSpaceTag> >;

  virtual ~AbstractLocalOperatorInterface() = default;

  virtual void alpha_volume(const EG&, const LFSU&, const X&, const LFSV&, R&) const
  {}

  virtual void alpha_boundary(const IG&, const LFSU&, const X&, const LFSV&, R&) const
  {}
};

#endif
