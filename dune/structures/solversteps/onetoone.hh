#ifndef DUNE_STRUCTURES_SOLVERSTEPS_ONETOONE_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_ONETOONE_HH

#include<dune/structures/onetoone.hh>
#include<dune/structures/solversteps/base.hh>

#include<memory>


template<int dim, typename... V>
class OneToOneMappingCheckerImpl
 : public TransitionSolverStepBase<V...>
{
  public:
  using Base = TransitionSolverStepBase<V...>;

  virtual ~OneToOneMappingCheckerImpl() {}

  void apply(std::shared_ptr<typename Base::Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    std::cout << "One-to-one checker not implemented for dimension " << dim << std::endl;
  }
};


template<typename... V>
class OneToOneMappingCheckerImpl<3, V...>
 : public TransitionSolverStepBase<V...>
{
  public:
  using Base = TransitionSolverStepBase<V...>;

  virtual ~OneToOneMappingCheckerImpl() {}

  void apply(std::shared_ptr<typename Base::Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    auto& gfs = vector->gridFunctionSpace();
    auto es = gfs.entitySet();

    Dune::PDELab::VectorDiscreteGridFunction vdgf(gfs, *vector);

    std::cout << "Checking one-to-one property of displacement field... ";
    std::cout << (is_onetoone(vdgf) ? "Success!" : "Failure") << std::endl;
  }
};


template<typename... V>
using OneToOneMappingChecker = OneToOneMappingCheckerImpl<TransitionSolverStepBase<V...>::dim, V...>;

#endif
