#ifndef DUNE_BLOCKLAB_OPERATORS_VIRTUALINTERFACE_HH
#define DUNE_BLOCKLAB_OPERATORS_VIRTUALINTERFACE_HH

/** An abstract base class for a local operator.
 *
 * It fixes all the duck-typed arguments in the process.
 * TODO: * Currently we are fixing numerical jacobians, which is not desirable.
 */

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/gridfunctionspace/localvector.hh>
#include<dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include<dune/pdelab/gridoperator/common/localmatrix.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/numericaljacobian.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/idefault.hh>

#include<memory>


namespace Dune::BlockLab {

  /** An abstract base class for a PDELab LocalOperator.
   *
   * The PDELab interface for a local operator uses templated signatures
   * for its methods which rely on duck-typing. This base class however
   * fixes all of the template parameters at the class level and therefore
   * allows to mark the methods virtual. This prohibits some advanced use
   * cases of PDELab, but it allows to select an operator at runtime!
   */
  template<typename GFSU, typename GFSV=GFSU>
  class AbstractLocalOperatorInterface
    : public Dune::PDELab::LocalOperatorDefaultFlags
    , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
    , public Dune::PDELab::NumericalJacobianVolume<AbstractLocalOperatorInterface<GFSU, GFSV>>
    , public Dune::PDELab::NumericalJacobianBoundary<AbstractLocalOperatorInterface<GFSU, GFSV>>
    , public Dune::PDELab::NumericalJacobianSkeleton<AbstractLocalOperatorInterface<GFSU, GFSV>>
    , public Dune::PDELab::FullVolumePattern
    , public Dune::PDELab::FullBoundaryPattern
    , public Dune::PDELab::FullSkeletonPattern
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

    enum { doAlphaBoundary = true };
    enum { doAlphaVolume = true };
    enum { doAlphaSkeleton = true };

    enum { doPatternVolume = true };
    enum { doPatternBoundary = true };
    enum { doPatternSkeleton = true };

    virtual void alpha_volume(const EG&, const LFSU&, const X&, const LFSV&, R&) const
    {}

    virtual void alpha_boundary(const IG&, const LFSU&, const X&, const LFSV&, R&) const
    {}

    virtual void alpha_skeleton(const IG&, const LFSU&, const X&, const LFSV&, const LFSU&, const X&, const LFSV&, R&, R&) const
    {}
  };

  /** A wrapper class that mixes the abstract interface into a LocalOperator class
   *
   * This can be wrapped around an existing local operator in order to make it
   * accessible through the virtual interface.
   */
  template<typename LocalOperator, typename GFSU, typename GFSV=GFSU>
  class VirtualizedLocalOperator
    : public AbstractLocalOperatorInterface<GFSU, GFSV>
  {
    using Base = AbstractLocalOperatorInterface<GFSU, GFSV>;
    using EG = typename Base::EG;
    using IG = typename Base::IG;
    using LFSU = typename Base::LFSU;
    using LFSV = typename Base::LFSV;
    using X = typename Base::X;
    using R = typename Base::R;

    public:
    VirtualizedLocalOperator(std::shared_ptr<LocalOperator> lop)
      : lop(lop)
    {}

    virtual ~VirtualizedLocalOperator() = default;

    virtual void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const override
    {
      if constexpr (LocalOperator::doAlphaVolume)
        lop->alpha_volume(eg, lfsu, x, lfsv, r);
    }

    virtual void alpha_boundary(const IG& ig, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const override
    {
      if constexpr (LocalOperator::doAlphaBoundary)
        lop->alpha_boundary(ig, lfsu, x, lfsv, r);
    }

    virtual void alpha_skeleton(const IG& ig, const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s, const LFSU& lfsu_n,
				              const X& x_n, const LFSV& lfsv_n, R& r_s, R& r_n) const override
    {
      if constexpr (LocalOperator::doAlphaSkeleton)
        lop->alpha_skeleton(ig, lfsu_s, x_s, lfsv_s, lfsu_n, x_n, lfsv_n, r_s, r_n);
    }

    private:
    std::shared_ptr<LocalOperator> lop;
  };

} // namespace Dune::BlockLab

#endif
