#ifndef DUNE_STRUCTURES_OPERATORINTERFACE_HH
#define DUNE_STRUCTURES_OPERATORINTERFACE_HH

/** An abstract base class for a local operator.
 *
 * It fixes all the duck-typed arguments in the process.
 */

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/gridfunctionspace/localvector.hh>
#include<dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include<dune/pdelab/gridoperator/common/localmatrix.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/structures/material.hh>


template<typename GFSU, typename GFSV>
class AbstractLocalOperatorInterface
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , public Dune::PDELab::NumericalJacobianVolume<AbstractLocalOperatorInterface<GFSU, GFSV>>
  , public Dune::PDELab::NumericalJacobianBoundary<AbstractLocalOperatorInterface<GFSU, GFSV>>
  , public Dune::PDELab::FullVolumePattern
  , public Dune::PDELab::FullBoundaryPattern
{
  public:
  using EG = Dune::PDELab::ElementGeometry<typename GFSU::Traits::GridView::template Codim<0>::Entity>;
  using IG = Dune::PDELab::IntersectionGeometry<typename GFSU::Traits::GridView::Intersection>;
  using LFSU = Dune::PDELab::LocalFunctionSpace<GFSU>;
  using LFSV = Dune::PDELab::LocalFunctionSpace<GFSV>;
  using X = Dune::PDELab::LocalVector<double, Dune::PDELab::TrialSpaceTag>;
  using R = Dune::PDELab::WeightedVectorAccumulationView<Dune::PDELab::LocalVector<double, Dune::PDELab::TestSpaceTag> >;
  using M = Dune::PDELab::WeightedMatrixAccumulationView<Dune::PDELab::LocalMatrix<double>>;

  virtual ~AbstractLocalOperatorInterface() = default;

  // These flags need to be put into the base class because they are otherwise sliced off.
  // Quite unfortunate (and preventing this base class from going into PDELab)
  enum { doAlphaBoundary = true };
  enum { doAlphaVolume = true };

  enum { doPatternVolume = true };
  enum { doPatternBoundary = true };

  virtual void alpha_volume(const EG&, const LFSU&, const X&, const LFSV&, R&) const
  {}

  virtual void alpha_boundary(const IG&, const LFSU&, const X&, const LFSV&, R&) const
  {}

  // I am not really happy to have this in here, but dont care too much
  virtual void setMaterial(std::shared_ptr<ElasticMaterialBase<typename GFSU::Traits::EntitySet, double>>)
  {}
};

#endif
