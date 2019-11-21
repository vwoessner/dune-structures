#ifndef DUNE_STRUCTURES_SOLVERSTEPS_TRANSFORMATION_HH
#define DUNE_STRUCTURES_SOLVERSTEPS_TRANSFORMATION_HH


#include<dune/pdelab.hh>
#include<dune/structures/solversteps/base.hh>
#include<dune/structures/solversteps/parametrizedlambda.hh>

#include<functional>


template<typename Vector>
class TransformationGridFunction
  : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<typename Vector::GridFunctionSpace::Traits::EntitySet,
                                       typename Vector::field_type,
                                       Dune::TypeTree::StaticDegree<typename Vector::GridFunctionSpace>::value,
                                       Dune::FieldVector<typename Vector::field_type, Dune::TypeTree::StaticDegree<typename Vector::GridFunctionSpace>::value> >,
                                       TransformationGridFunction<Vector>
      >
{
  public:
  using Traits = Dune::PDELab::GridFunctionTraits<typename Vector::GridFunctionSpace::Traits::EntitySet,
      typename Vector::field_type,
      Dune::TypeTree::StaticDegree<typename Vector::GridFunctionSpace>::value,
      Dune::FieldVector<typename Vector::field_type, Dune::TypeTree::StaticDegree<typename Vector::GridFunctionSpace>::value> >;

  using GridFunctionSpace = typename Vector::GridFunctionSpace;
  using GridFunction =  typename Dune::PDELab::VectorDiscreteGridFunction<GridFunctionSpace, Vector>;


  TransformationGridFunction(std::function<typename Traits::RangeType(typename Traits::DomainType, typename Traits::RangeType)> func,
                             Vector& vector)
    : func(func), gridfunction(vector.gridFunctionSpace(), vector)
  {}

  inline void evaluate(const typename Traits::ElementType& e,
                       const typename Traits::DomainType& x,
                       typename Traits::RangeType& y) const
  {
    typename Traits::RangeType eval;
    gridfunction.evaluate(e, x, eval);
    auto xg = e.geometry().global(x);
    y = func(eval, xg);
  }

  private:
  GridFunction gridfunction;
  std::function<typename Traits::RangeType(typename Traits::DomainType, typename Traits::RangeType)> func;
};


template<typename Vector>
class TransformationTransitionStep : public TransitionSolverStepBase<Vector>
{
  public:
  using Base = TransitionSolverStepBase<Vector>;
  using GridFunction = TransformationGridFunction<Vector>;

  using FunctionSignature = Dune::FieldVector<double, 3>(Dune::FieldVector<double, 3>, Dune::FieldVector<double, 3>);

  template<typename FUNC>
  TransformationTransitionStep(const FUNC& func)
    : func(func)
  {}

  virtual ~TransformationTransitionStep() {}

  virtual void apply(Vector& vector, typename Base::ConstraintsContainer&) override
  {
    std::cout << "Transforming solution!" << std::endl;
    TransformationGridFunction<Vector> trafo(func, vector);
    Dune::PDELab::interpolate(trafo, vector.gridFunctionSpace(), vector);
  }

  private:
  std::function<FunctionSignature> func;
};


template<typename Vector, typename... P>
using ParametrizedTransformationTransitionStep = ParametrizedLambdaVariationTransitionStepBase<TransformationTransitionStep<Vector>, P...>;

#endif
