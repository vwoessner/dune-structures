#ifndef DUNE_STRUCTURES_VON_MISES_HH
#define DUNE_STRUCTURES_VON_MISES_HH

#include<dune/pdelab.hh>
#include<dune/structures/material.hh>

#include<memory>

/** Calculate von Mises Stress */
template<typename Container>
class VonMisesStressGridFunction
		: public Dune::PDELab::GridFunctionBase<
      Dune::PDELab::GridFunctionTraits<typename Container::GridFunctionSpace::Traits::GridViewType,
                                       typename Container::field_type,
                                       1,
                                       Dune::FieldVector<typename Container::field_type,1> >,
      VonMisesStressGridFunction<Container>
			>
{
  private:
	using DisplacementFunction = Dune::PDELab::VectorDiscreteGridFunction<typename Container::GridFunctionSpace, Container>;
	using DisplacementFunctionGradient = Dune::PDELab::VectorDiscreteGridFunctionGradient<typename Container::GridFunctionSpace, Container>;
	using Material = ElasticMaterialBase<typename Container::GridFunctionSpace::Traits::GridViewType, typename Container::field_type>;

  public:
  using Traits = Dune::PDELab::GridFunctionTraits<
  		typename Container::GridFunctionSpace::Traits::GridViewType,
      typename Container::field_type,
      1,
      Dune::FieldVector<typename Container::field_type,1>
  >;

	VonMisesStressGridFunction(const Container& cont, std::shared_ptr<Material> material)
    : displacement_grad_f(cont.gridFunctionSpace(), cont),
			material(material)
  {}

  inline void evaluate(const typename Traits::ElementType& e,
  		                 const typename Traits::DomainType& x,
											 typename Traits::RangeType& y) const
  {
  	// Evaluate the displacement and its gradient
  	typename DisplacementFunctionGradient::Traits::RangeType displacement_grad;
  	displacement_grad_f.evaluate(e, x, displacement_grad);

  	// Get the Lame coefficients for this problem
  	auto lame1 = material->first_lame(e, x);
  	auto lame2 = material->second_lame(e, x);

  	typename Traits::RangeType sum = 0.0;

  	// Offdiagonal part of the deviatoric stress tensor
  	for (int i=0; i<2; ++i)
  		for (int j=0; j<i; ++j)
  		{
  		  auto entry = lame2 * (displacement_grad[i][j] + displacement_grad[j][i]);
  			sum += 2.0 * entry * entry;
  		}

  	// Calculate divergence
  	typename Traits::RangeType div = 0.0;
  	for (int i=0; i<2; ++i)
  		div += displacement_grad[i][i];

  	std::array<typename Traits::RangeType, 2> sig;
  	for(int i=0; i<2; ++i)
  	  sig[i] = lame1 * div + 2.0 * lame2 * displacement_grad[i][i];

  	// Diagonal part of the deviatoric stress tensor
  	for (int i=0; i<2; ++i)
  		for (int j=0; j<2; ++j)
  		{
  			auto entry = (i == j ? 2./3. : -1./3.) * sig[j];
  			sum += entry * entry;
  		}

  	using std::sqrt;
  	y = std::sqrt(1.5 * sum);
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return displacement_grad_f.getGridView();
  }

  private:
  DisplacementFunctionGradient displacement_grad_f;
  std::shared_ptr<Material> material;
};

#endif
