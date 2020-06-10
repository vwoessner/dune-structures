#ifndef DUNE_BLOCKLAB_GRIDS_STRUCTURED_HH
#define DUNE_BLOCKLAB_GRIDS_STRUCTURED_HH

#include<dune/blocklab/grids/capabilities.hh>
#include<dune/common/filledarray.hh>
#include<dune/common/fvector.hh>
#include<dune/common/parametertree.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/pdelab/common/partitionviewentityset.hh>

#include<memory>
#include<tuple>


namespace Dune::BlockLab {

  namespace impl {

    template<typename GridImpl, bool simplex>
    class StructuredGridProviderImpl
    {
      public:
      using Grid = GridImpl;
      using Parameter = std::tuple<>;

      static constexpr int dim = Grid::dimension;

      StructuredGridProviderImpl(const Dune::ParameterTree& config)
	: config(config)
      {}

      std::shared_ptr<Grid> createGrid()
      {
	if(!grid)
	{
	  auto ll = config.get<Dune::FieldVector<double, dim>>("lowerleft", Dune::FieldVector<double, dim>(0.0));
	  auto ur = config.get<Dune::FieldVector<double, dim>>("upperright", Dune::FieldVector<double, dim>(1.0));
	  auto N = config.get<std::array<unsigned int, dim>>("N", Dune::filledArray<dim, unsigned int>(10));

	  if constexpr (simplex)
	    grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(ll, ur, N);
	  else
	    grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(ll, ur, N);
	  grid->loadBalance();
	  grid->globalRefine(config.get<int>("refinement", 0));
	}

	return grid;
      }

      std::shared_ptr<Parameter> createParameters()
      {
	return std::make_shared<Parameter>();
      }

      private:
      std::shared_ptr<Grid> grid;
      Dune::ParameterTree config;
    };

  } // namespace impl


  template<typename Grid>
  class StructuredSimplexGridProvider
    : public impl::StructuredGridProviderImpl<Grid, true>
  {
    using impl::StructuredGridProviderImpl<Grid, true>::StructuredGridProviderImpl;
  };

  template<typename Grid>
  class StructuredCubeGridProvider
    : public impl::StructuredGridProviderImpl<Grid, false>
  {
    using impl::StructuredGridProviderImpl<Grid, false>::StructuredGridProviderImpl;
  };

  namespace Capabilities {

    template<typename Grid>
    struct HasSimplices<StructuredSimplexGridProvider<Grid>>
    {
      static constexpr bool value = true;
    };

    template<typename Grid>
    struct HasCubes<StructuredCubeGridProvider<Grid>>
    {
      static constexpr bool value = true;
    };

    template<typename Grid>
    struct HasSimplices<StructuredCubeGridProvider<Grid>>
    {
      static constexpr bool value = isAdaptive<StructuredCubeGridProvider<Grid>>();
    };

  } // namespace Capabilities

} // namespace Dune::BlockLab

#endif
