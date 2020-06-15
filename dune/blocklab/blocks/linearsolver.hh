#ifndef DUNE_BLOCKLAB_BLOCKS_LINEARSOLVER_HH
#define DUNE_BLOCKLAB_BLOCKS_LINEARSOLVER_HH

/* A solver block for a linear solver.
 * This currently unconditionally uses a direct solver, which is
 * completely unacceptable. However, configurable ISTL solvers are
 * currently implemented upstream and I want to transition directly
 * to that work instead of building another thing.
 */

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/operators/virtualinterface.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/common/parametertree.hh>

#include<memory>

namespace Dune::BlockLab {

  template<typename P, typename V, std::size_t i>
  class LinearSolverBlock
     : public BlockBase<P, V, i>
  {
    public:
    using Traits = BlockTraits<P, V, i>;

    using VirtLocalOperator = AbstractLocalOperatorInterface<typename Traits::GridFunctionSpace>;
    using GridOperator = Dune::PDELab::GridOperator<typename Traits::GridFunctionSpace,
						    typename Traits::GridFunctionSpace,
						    VirtLocalOperator,
						    Dune::PDELab::ISTL::BCRSMatrixBackend<>,
						    typename Traits::ctype,
						    typename Traits::Range,
						    typename Traits::Range,
						    typename Traits::ConstraintsContainer,
						    typename Traits::ConstraintsContainer>;

    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_UMFPack;
    using StationaryLinearProblemSolver = Dune::PDELab::StationaryLinearProblemSolver<GridOperator, LinearSolver, typename Traits::Vector>;

    template<typename Context>
    LinearSolverBlock(Context&, const Dune::ParameterTree& config)
      : operatorstr(config.get<std::string>("operator"))
    {}

    virtual ~LinearSolverBlock() = default;

    virtual void update_parameter(std::string name, typename Traits::Parameter param) override
    {
      if (name == operatorstr)
      {
	// The operator changed - rebuild the linear solver object!
	auto vector = this->solver->template getVector<i>();
	auto cc = this->solver->template getConstraintsContainer<i>();
	auto gfs = vector->gridFunctionSpaceStorage();
	auto localoperator = this->solver->template param<std::shared_ptr<VirtLocalOperator>>(operatorstr);
	Dune::PDELab::ISTL::BCRSMatrixBackend<> mb(21);
	gridoperator = std::make_shared<GridOperator>(*gfs, *cc, *gfs, *cc, *localoperator, mb);
	linearsolver = std::make_shared<LinearSolver>(0);
	slp = std::make_shared<StationaryLinearProblemSolver>(*gridoperator, *linearsolver, *vector, 1e-12);
      }
    }

    virtual void apply() override
    {
      std::cout << "Applying linear Solver!" << std::endl;
      slp->apply();
    }

    protected:
    std::string operatorstr;
    std::shared_ptr<LinearSolver> linearsolver;
    std::shared_ptr<GridOperator> gridoperator;
    std::shared_ptr<StationaryLinearProblemSolver> slp;
  };

} // namespace Dune::BlockLab

#endif
