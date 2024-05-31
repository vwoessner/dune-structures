#ifndef DUNE_STRUCTURES_VON_MISES_HH
#define DUNE_STRUCTURES_VON_MISES_HH

#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/structures/material.hh>

#include <array>
#include <memory>

/// Calculate von Mises Stress
template<typename Container, int dim>
class VonMisesStressGridFunction
  : public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<
        typename Container::GridFunctionSpace::Traits::EntitySet,
        typename Container::field_type,
        1,
        Dune::FieldVector<typename Container::field_type, 1>>,
      VonMisesStressGridFunction<Container, dim>>
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
    typename Container::GridFunctionSpace::Traits::GridView,
    typename Container::field_type,
    1,
    Dune::FieldVector<typename Container::field_type, 1>>;

  VonMisesStressGridFunction(const Container& cont,
                             std::shared_ptr<Material> material)
    : displacement_grad_f(cont.gridFunctionSpace(), cont)
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

    using RF = typename Traits::RangeType;

    // Calculate the Cauchy stress tensor
    Dune::FieldMatrix<RF, dim, dim> stress;
    for (int i = 0; i < dim; ++i)
    {
      for (int j = 0; j <= i; ++j)
      {
        const RF val =
          lame2 * (displacement_grad[i][j] + displacement_grad[j][i]);
        stress[i][j] = val;
        stress[j][i] = val;
      }
    }
    RF trace = 0.0;
    for (int i = 0; i < dim; ++i)
    {
      trace += displacement_grad[i][i];
    }
    for (int i = 0; i < dim; ++i)
    {
      stress[i][i] += lame1 * trace;
    }

    // Calculate von Mises stress
    RF sum = 0.0;
    // Diagonal
    for (int i = 0; i < dim; ++i)
    {
      sum += stress[i][i] * stress[i][i];
    }
    // Off-diagonal
	if (dim == 2)
	{
		sum += 3 * stress[0][1] * stress[0][1] - stress[0][0] * stress[1][1];
	}
	else if (dim == 3)
	{
		for (int i = 0; i < dim; ++i)
		{
		  int j = i + 1;
		  if (i == dim - 1)
			j = 0;

		  sum += 3 * stress[i][j] * stress[i][j] - stress[i][i] * stress[j][j];
		}
	}
	else
	{
		std::cout << "Von mises stress makes no sense in the given dimension";
	}
    y = std::sqrt(sum);
  }

  inline const typename Traits::GridViewType& getGridView() const
  {
    return displacement_grad_f.getGridView();
  }

private:
  DisplacementFunctionGradient displacement_grad_f;
  std::shared_ptr<Material> material;
};

#endif
