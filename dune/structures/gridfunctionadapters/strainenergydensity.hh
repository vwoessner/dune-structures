#ifndef DUNE_STRUCTURES_STRAINENERGYDENSITY_HH
#define DUNE_STRUCTURES_STRAINENERGYDENSITY_HH

#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/structures/material.hh>

/// Calculate the strain energy density
template<class Vector>
class StrainEnergyDensityFunction
  : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<
        typename Vector::GridFunctionSpace::Traits::EntitySet,
        typename Vector::field_type,
        1,
        Dune::FieldVector<typename Vector::field_type, 1>>,
      StrainEnergyDensityFunction<Vector>>
{
private:
  using DisplacementFunction =
    Dune::PDELab::VectorDiscreteGridFunction<typename Vector::GridFunctionSpace,
                                             Vector>;
  using DisplacementFunctionGradient =
    Dune::PDELab::VectorDiscreteGridFunctionGradient<
      typename Vector::GridFunctionSpace,
      Vector>;
  using Material =
    ElasticMaterialBase<typename Vector::GridFunctionSpace::Traits::EntitySet,
                        typename Vector::field_type>;

public:
  using Traits = Dune::PDELab::GridFunctionTraits<
    typename Vector::GridFunctionSpace::Traits::GridView,
    typename Vector::field_type,
    1,
    Dune::FieldVector<typename Vector::field_type, 1>>;

  StrainEnergyDensityFunction(const Vector& vec,
                              std::shared_ptr<Material> material)
    : displacement_grad_f(vec.gridFunctionSpace(), vec)
    , material(material)
  {
  }

  inline void evaluate(const typename Traits::ElementType& e,
                       const typename Traits::DomainType& x,
                       typename Traits::RangeType& y) const
  {
    // Evaluate the displacement and its gradient
    typename DisplacementFunctionGradient::Traits::RangeType displacement_grad;
    displacement_grad_f.evaluate(e, x, displacement_grad);

    // Get the Lame coefficients for this problem
    auto lame1 = material->parameter(e, x, 0);
    auto lame2 = material->parameter(e, x, 1);

    // Calculate strain tensor
    using RF = typename Traits::RangeType;
    constexpr size_t dim =
      Vector::GridFunctionSpace::Traits::GridView::dimension;
    Dune::FieldMatrix<RF, dim, dim> strain;
    for (size_t i = 0; i < dim; ++i)
    {
      for (size_t j = 0; j <= i; ++j)
      {
        const RF val =
          0.5 * (displacement_grad[i][j] + displacement_grad[j][i]);
        strain[i][j] = val;
        strain[j][i] = val;
      }
    }

    // Calculate strain energy density
    RF trace(0.0);
    for (size_t i = 0; i < dim; ++i)
    {
      trace += strain[i][i];
    }

    RF entries_squared(0.0);
    for (size_t i = 0; i < dim; ++i)
    {
      for (size_t j = 0; j < dim; ++j)
      {
        entries_squared += strain[i][j] * strain[i][j];
      }
    }

    y = 0.5 * lame1 * trace * trace + lame2 * entries_squared;
  }

  inline const typename Traits::GridViewType& getGridView() const
  {
    return displacement_grad_f.getGridView();
  }

private:
  DisplacementFunctionGradient displacement_grad_f;
  std::shared_ptr<Material> material;
};
#endif // DUNE_STRUCTURES_STRAINENERGYDENSITY_HH
