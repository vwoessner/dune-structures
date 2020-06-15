#ifndef DUNE_BLOCKLAB_BLOCKS_CONSTRAINTS_HH
#define DUNE_BLOCKLAB_BLOCKS_CONSTRAINTS_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/construction/callabletree.hh>
#include<dune/blocklab/construction/muparser.hh>
#include<dune/common/parametertree.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/typetree/utility.hh>

#include<array>
#include<memory>


namespace Dune::BlockLab {

  template<typename P, typename V, std::size_t i>
  class ConstraintsBlock
    : public BlockBase<P, V, i>
  {
    public:
    using Traits = typename BlockBase<P, V, i>::Traits;
    using FunctionSignature = bool(typename Traits::GridView::Intersection, typename Traits::GridView::Intersection::Geometry::LocalCoordinate);

    template<typename Context>
    ConstraintsBlock(Context& ctx, const Dune::ParameterTree& config)
      : ConstraintsBlock(muparser_callable_array<Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount, FunctionSignature>(
                         config.get<std::string>("functions"),
                         ctx.getSolver())
			 )
    {}

    ConstraintsBlock(const std::array<std::function<FunctionSignature>,
                                      Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount
                                      >& funcs)
      : funcs(funcs)
    {}

    virtual ~ConstraintsBlock() {}

    virtual void apply() override
    {
      auto constraintscontainer = this->solver->template getConstraintsContainer<i>();
      auto gfs = this->solver->template getVector<i>()->gridFunctionSpaceStorage();
      auto bctype = makeBoundaryConditionTreeFromCallables(*gfs, funcs);

      Dune::PDELab::constraints(bctype, *gfs, *constraintscontainer);
      std::cout << "Assembled constraints - " << constraintscontainer->size() << " of " << gfs->size() << " dofs constrained!" << std::endl;
    }

    private:
    // Store the lambdas
    std::array<std::function<FunctionSignature>,
               Dune::TypeTree::TreeInfo<typename Traits::GridFunctionSpace>::leafCount
               > funcs;
  };

} // namespace Dune::BlockLab

#endif
