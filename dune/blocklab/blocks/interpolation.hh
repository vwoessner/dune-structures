#ifndef DUNE_BLOCKLAB_BLOCKS_INTERPOLATION_HH
#define DUNE_BLOCKLAB_BLOCKS_INTERPOLATION_HH

/** A solver block that interpolates a function given through
 *  an array (one entry per component of the function space tree)
 *  into a DoF vector. Context-based construction uses the
 *  MuParser library to parse functions from the inifile.
 */
#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/construction/callabletree.hh>
#include<dune/blocklab/construction/muparser.hh>
#include<dune/common/parametertree.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/typetree/utility.hh>

#include<array>
#include<functional>
#include<iostream>
#include<memory>


namespace Dune::BlockLab {

  template<typename P, typename V, std::size_t i>
  class InterpolationBlock
    : public BlockBase<P, V, i>
  {
    public:
    using Traits = typename BlockBase<P, V, i>::Traits;
    using FunctionSignature = typename Traits::Range(typename Traits::Entity, typename Traits::GlobalCoordinate);

    /* Construct this block from the construction context */
    template<typename Context>
    InterpolationBlock(Context& ctx, const Dune::ParameterTree& config)
      : InterpolationBlock(muparser_callable_array<Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount, FunctionSignature>(
	                    config.get<std::string>("functions"),
                            ctx.getSolver())
			   )
    {}

    /** Construct this block from an array of callables */
    InterpolationBlock(const std::array<std::function<FunctionSignature>,
                                        Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount
                                        >& funcs)
      : funcs(funcs)
    {}

    virtual ~InterpolationBlock() = default;

    virtual void apply() override
    {
      auto vector = this->solver->template getVector<i>();
      auto gfs = vector->gridFunctionSpaceStorage();
      auto gf = Dune::BlockLab::makeGridFunctionTreeFromCallables(*gfs, funcs);

      std::cout << "Interpolating into solution vector" << std::endl;
      Dune::PDELab::interpolate(gf, *gfs, *vector);
    }

    private:
    // Store the lambdas
    std::array<std::function<FunctionSignature>,
	       Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount
	       > funcs;
  };

} // namespace Dune::BlockLab

#endif
