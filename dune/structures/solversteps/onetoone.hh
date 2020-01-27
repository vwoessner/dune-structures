#ifndef DUNE_STRUCTURES_SOLVERSTEPS_ONETOONE_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_ONETOONE_HH

#include<dune/structures/onetoone.hh>
#include<dune/structures/solversteps/base.hh>

#include<memory>


template<typename Vector, int dim>
class OneToOneMappingCheckerImpl
 : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  virtual ~OneToOneMappingCheckerImpl() {}

  void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    std::cout << "One-to-one checker not implemented for dimension " << dim << std::endl;
  }
};


template<typename Vector>
class OneToOneMappingCheckerImpl<Vector, 3>
 : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;

  virtual ~OneToOneMappingCheckerImpl() {}

  void apply(std::shared_ptr<Vector> vector, std::shared_ptr<typename Base::ConstraintsContainer> cc) override
  {
    auto& gfs = vector->gridFunctionSpace();
    auto es = gfs.entitySet();

    Dune::PDELab::VectorDiscreteGridFunction vdgf(gfs, *vector);

    std::cout << "Checking one-to-one property of displacement field... ";
    std::cout << (is_onetoone(vdgf) ? "Success!" : "Failure") << std::endl;
  }
};


template<typename Vector>
using OneToOneMappingChecker = OneToOneMappingCheckerImpl<Vector, TransitionSolverStepBase<Vector>::dim>;

#endif
