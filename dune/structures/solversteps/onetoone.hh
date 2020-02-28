#ifndef DUNE_STRUCTURES_SOLVERSTEPS_ONETOONE_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_ONETOONE_HH

#include<dune/structures/onetoone.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/traits.hh>

#include<memory>


template<int dim, typename... V>
class OneToOneMappingCheckerImpl
 : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  virtual ~OneToOneMappingCheckerImpl() {}

  void apply(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc) override
  {
    std::cout << "One-to-one checker not implemented for dimension " << dim << std::endl;
  }
};


template<typename... V>
class OneToOneMappingCheckerImpl<3, V...>
 : public TransitionSolverStepBase<V...>
{
  public:
  using Traits = SimpleStepTraits<V...>;

  virtual ~OneToOneMappingCheckerImpl() {}

  void apply(std::shared_ptr<typename Traits::Vector> vector, std::shared_ptr<typename Traits::ConstraintsContainer> cc) override
  {
    auto& gfs = vector->gridFunctionSpace();
    auto es = gfs.entitySet();

    Dune::PDELab::VectorDiscreteGridFunction vdgf(gfs, *vector);

    std::cout << "Checking one-to-one property of displacement field... ";
    std::cout << (is_onetoone(vdgf) ? "Success!" : "Failure") << std::endl;
  }
};


template<typename... V>
using OneToOneMappingChecker = OneToOneMappingCheckerImpl<SimpleStepTraits<V...>::dim, V...>;

#endif
