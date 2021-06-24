#ifndef DUNE_STRUCTURES_STRESS_HH
#define DUNE_STRUCTURES_STRESS_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fmatrixev.hh>
#include <dune/common/fvector.hh>

#include <dune/pdelab.hh>
#include <dune/structures/material.hh>

#include <memory>

/// Calculate the stress eigenvectors and evaluate a selected one
/**
 *  The evaluate() method calculates the local stress matrix and its
 *  eigenvectors and eigenvalues. They are ordered by size. The eigenvector
 *  with the index selected via set_index() is returned.
 */
template<typename Container, int dim>
class StressEVGridFunction
  : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<
        typename Container::GridFunctionSpace::Traits::EntitySet,
        typename Container::field_type,
        dim,
        Dune::FieldVector<typename Container::field_type, dim>>,
      StressEVGridFunction<Container, dim>>
{
private:
  using DisplacementFunction = Dune::PDELab::VectorDiscreteGridFunction<
    typename Container::GridFunctionSpace,
    Container>;
  using DisplacementFunctionGradient =
    Dune::PDELab::VectorDiscreteGridFunctionGradient<
      typename Container::GridFunctionSpace,
      Container>;
  using Material = ElasticMaterialBase<
    typename Container::GridFunctionSpace::Traits::EntitySet,
    typename Container::field_type>;

public:
  using Traits = Dune::PDELab::GridFunctionTraits<
    typename Container::GridFunctionSpace::Traits::EntitySet,
    typename Container::field_type,
    dim,
    Dune::FieldVector<typename Container::field_type, dim>>;

public:
  StressEVGridFunction(const Container& cont,
                       const std::shared_ptr<Material> material)
    : displacement_grad_f(cont.gridFunctionSpace(), cont)
    , _material(material)
  {
  }

  inline void evaluate(const typename Traits::ElementType& e,
                       const typename Traits::DomainType& x,
                       typename Traits::RangeType& y) const
  {
    if (not _material)
      DUNE_THROW(Dune::Exception, "Material pointer not valid!");

    // Evaluate the displacement and its gradient
    typename DisplacementFunctionGradient::Traits::RangeType displacement_grad;
    displacement_grad_f.evaluate(e, x, displacement_grad);

    // Evaluate the strain tensor
    using RF = typename Traits::RangeFieldType;
    Dune::FieldMatrix<RF, dim, dim> strain;
    RF strain_trace{ 0.0 };
    for (size_t i = 0; i < dim; ++i)
    {
      for (size_t j = 0; j <= i; ++j)
      {
        strain[i][j] =
          0.5 * (displacement_grad[i][j] + displacement_grad[j][i]);
        if (i == j)
          strain_trace += strain[i][i];
        else
          strain[j][i] = strain[i][j];
      }
    }

    // Get the Lame coefficients for this problem
    const auto lame1 = _material->parameter(e, x, 0); // \lambda
    const auto lame2 = _material->parameter(e, x, 1); // \mu

    // Evaluate stress tensor
    Dune::FieldMatrix<RF, dim, dim> stress;
    for (size_t i = 0; i < dim; ++i)
    {
      for (size_t j = 0; j <= i; ++j)
      {
        stress[i][j] = 2 * lame2 * strain[i][j];
        if (i == j)
        {
          stress[i][i] += lame1 * strain_trace;
        }
        else
        {
          stress[j][i] = stress[i][j];
        }
      }
    }

    // Compute Eigenvectors and Eigenvalues
    Dune::FieldMatrix<RF, dim, dim> eigenvectors;
    Dune::FieldVector<RF, dim> eigenvalues;
    Dune::FMatrixHelp::eigenValuesVectors(stress, eigenvalues, eigenvectors);

    // Output selected eigenvector
    y = eigenvectors[_vector_idx] * eigenvalues[_vector_idx];
  }

  inline const typename Traits::GridViewType& getGridView() const
  {
    return displacement_grad_f.getGridView();
  }

  /// Select the index of the eigenvector to evaluate
  void set_index(const size_t index) { _vector_idx = index; }

private:
  DisplacementFunctionGradient displacement_grad_f;
  std::shared_ptr<Material> _material;
  /// Index of the eigenvector to plot. Defaults to largest eigenvalue.
  size_t _vector_idx = dim - 1;
};

#endif // DUNE_STRUCTURES_STRESS_HH
