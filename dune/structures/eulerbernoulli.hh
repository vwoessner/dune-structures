#ifndef DUNE_STRUCTURES_EULERBERNOULLI_HH
#define DUNE_STRUCTURES_EULERBERNOULLI_HH

/** Local operators to implement the beam-truss elements for
 *  Euler-Bernoulli beams used in
 *
 *  Hansbo P., Larson M. and Larsson K.
 *  Cut Finite Element Methods for Linear Elasticity Problems
 */

#include<dune/pdelab.hh>
#include<dune/pdelab/localoperator/sum.hh>
#include<dune/structures/elasticity.hh>
#include<dune/structures/material.hh>
#include<tuple>


/** An operator that implement a single Euler-Bernoulli beam
 *
 * For multiple beams, use a sum of operators. This may not be super
 * high performance, but things are already complicated enough!
 */
template<typename GFS>
class EulerBernoulli2DLocalOperator
  :  public Dune::PDELab::FullSkeletonPattern
  ,  public Dune::PDELab::FullVolumePattern
  ,  public Dune::PDELab::LocalOperatorDefaultFlags
  ,  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  ,  public Dune::PDELab::NumericalJacobianVolume<EulerBernoulli2DLocalOperator<GFS>>
  ,  public Dune::PDELab::NumericalJacobianSkeleton<EulerBernoulli2DLocalOperator<GFS>>
{
  public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  enum { doPatternSkeleton = true };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doAlphaSkeleton  = true };

  EulerBernoulli2DLocalOperator(const Dune::ParameterTree& params)
    : Dune::PDELab::NumericalJacobianVolume<EulerBernoulli2DLocalOperator<GFS>>(1e-9)
    , Dune::PDELab::NumericalJacobianSkeleton<EulerBernoulli2DLocalOperator<GFS>>(1e-9)
  {
    std::cout << "Constructed" << std::endl;
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    std::cout << "alpha_volume" << std::endl;
  }

  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton (const IG& ig,
                       const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                       const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                       R& r_s, R& r_n) const
  {
    std::cout << "alpha_skeleton" << std::endl;
  }
};


template<typename GFS, int dim>
class FibreReinforcedBulkOperator
  : public Dune::PDELab::InstationarySumLocalOperator<
             std::tuple<typename OperatorSwitch<GFS, dim>::Elasticity,
                        EulerBernoulli2DLocalOperator<GFS>>>
{
  public:
  using BulkOperator = typename OperatorSwitch<GFS, dim>::Elasticity;
  using FibreOperator = EulerBernoulli2DLocalOperator<GFS>;
  using BaseT = Dune::PDELab::InstationarySumLocalOperator<std::tuple<BulkOperator, FibreOperator>>;

  FibreReinforcedBulkOperator(std::shared_ptr<const GFS> gfs,
                              const Dune::ParameterTree& rootparams,
                              const Dune::ParameterTree& params,
                              std::shared_ptr<ElasticMaterialBase<typename GFS::Traits::EntitySet, double>> material)
    : bulkoperator(*gfs, *gfs, params, material)
    , fibreoperator(params)
  {
    this->template setSummand<0>(bulkoperator);
    this->template setSummand<1>(fibreoperator);
  }

  template<typename CGFS, typename CVEC>
  void setCoefficientForce(std::shared_ptr<CGFS> gfs, std::shared_ptr<CVEC> force)
  {
    bulkoperator.setCoefficientForce(gfs, force);
  }

  void setMaterial(std::shared_ptr<ElasticMaterialBase<typename GFS::Traits::EntitySet, double>> material)
  {
    bulkoperator.setMaterial(material);
  }

  template<typename CGFS, typename CVEC>
  void setCoefficientTraction(std::shared_ptr<CGFS> gfs, std::shared_ptr<CVEC> traction)
  {
    bulkoperator.setCoefficientTraction(gfs, traction);
  }

  private:
  typename OperatorSwitch<GFS, dim>::Elasticity bulkoperator;
  EulerBernoulli2DLocalOperator<GFS> fibreoperator;
};


#endif
