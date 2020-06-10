#ifndef DUNE_BLOCKLAB_BLOCKS_ERROR_HH
#define DUNE_BLOCKLAB_BLOCKS_ERROR_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/construction/callabletree.hh>
#include<dune/blocklab/construction/muparser.hh>
#include<dune/common/parametertree.hh>
#include<dune/pdelab/common/functionutilities.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionadapter.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/typetree/utility.hh>

#include<array>
#include<functional>
#include<memory>
#include<string>


namespace Dune::BlockLab {

  template<typename P, typename V, std::size_t i>
  class DiscretizationErrorBlock
    : public BlockBase<P, V, i>
  {
    public:
    using Traits = typename BlockBase<P, V, i>::Traits;
    using FunctionSignature = typename Traits::Range(typename Traits::Entity, typename Traits::GlobalCoordinate);

    /* Construct this block from the construction context */
    template<typename Context>
    DiscretizationErrorBlock(Context& ctx, const Dune::ParameterTree& config)
      : DiscretizationErrorBlock(muparser_callable_array<Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount, FunctionSignature>(
	                           config.get<std::string>("analytic"),
                                   ctx.getSolver()),
				 config.get<std::string>("norm", "l2"),
				 config.get<std::string>("name", "error")
			         )
    {}

    /** Construct this block from an array of callables and a norm string*/
    DiscretizationErrorBlock(const std::array<std::function<FunctionSignature>,
                                              Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount
                                              >& funcs,
                             std::string norm,
                             std::string name)
      : funcs(funcs)
      , norm(norm)
      , name(name)
    {}

    virtual ~DiscretizationErrorBlock() = default;

    virtual void setup() override
    {
      this->solver->template introduce_parameter<double>(name, 0.0);
    }

    virtual void apply() override
    {
      auto vector = this->solver->template getVector<i>();
      auto gfs = vector->gridFunctionSpaceStorage();

      // Make a grid function object for the analytic solution
      auto f_analytic = Dune::BlockLab::makeGridFunctionTreeFromCallables(*gfs, funcs);

      if (norm == "l2")
      {
	 Dune::PDELab::DiscreteGridFunction<typename Traits::GridFunctionSpace, typename Traits::Vector> dgf(gfs, vector);
	 Dune::PDELab::DifferenceSquaredAdapter diff(dgf, f_analytic);
	 // TODO Vector-valued stuff eternally broken
	 Dune::FieldVector<double, 1> sum = 0.0;
	 // TODO configure integration order
	 Dune::PDELab::integrateGridFunction(diff, sum, 4);
	 sum = dgf.getGridView().comm().sum(sum);

	 this->solver->update_parameter(name, std::sqrt(sum));
      }
      else
	DUNE_THROW(Dune::Exception, "Norm unknown in discretization error calculation!");
    }

    private:
    // Store the lambdas that describe the analytic solution
    std::array<std::function<FunctionSignature>,
	       Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount
	       > funcs;

    std::string norm;
    std::string name;
  };

} // namespace Dune::BlockLab

#endif
