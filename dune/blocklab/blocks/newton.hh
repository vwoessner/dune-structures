#ifndef DUNE_BLOCKLAB_BLOCKS_NEWTON_HH
#define DUNE_BLOCKLAB_BLOCKS_NEWTON_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/operators/virtualinterface.hh>
#include<dune/common/parametertree.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/solver/newton.hh>

#include<memory>

namespace Dune::BlockLab {

  template<typename P, typename V, std::size_t i>
  class NewtonSolverBlock
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
    using NewtonSolver = Dune::PDELab::NewtonMethod<GridOperator, LinearSolver>;

    template<typename Context>
    NewtonSolverBlock(Context&, const Dune::ParameterTree& config)
      : NewtonSolverBlock(config)
    {}

    NewtonSolverBlock(const Dune::ParameterTree& config)
      : config(config)
    {}

    virtual ~NewtonSolverBlock() {}

    virtual void update_parameter(std::string name, typename Traits::Parameter param) override
    {
      if (name == config.get<std::string>("operator"))
      {
        // The operator changed - rebuild the Newton solver object!
        auto vector = this->solver->template getVector<i>();
        auto cc = this->solver->template getConstraintsContainer<i>();
        auto gfs = vector->gridFunctionSpaceStorage();
        auto localoperator = this->solver->template param<std::shared_ptr<VirtLocalOperator>>(config.get<std::string>("operator"));
        Dune::PDELab::ISTL::BCRSMatrixBackend<> mb(21);
        gridoperator = std::make_shared<GridOperator>(*gfs, *cc, *gfs, *cc, *localoperator, mb);
        linearsolver = std::make_shared<LinearSolver>(0);
        std::cout << "gridoperator " << gridoperator.get() << std::endl;
        std::cout << "linearsolver " << linearsolver.get() << std::endl;
        newton = std::make_shared<NewtonSolver>(*gridoperator, *linearsolver, config);
        std::cout << "REBUILDING NEWTON" << std::endl;
        newton->setVerbosityLevel(2);
      }
    }

    virtual void apply() override
    {
      std::cout << "Applying Newton Solver!" << std::endl;
      newton->apply(*(this->solver->template getVector<i>()));
    }

    protected:
    Dune::ParameterTree config;
    std::shared_ptr<LinearSolver> linearsolver;
    std::shared_ptr<GridOperator> gridoperator;
    std::shared_ptr<NewtonSolver> newton;
  };

} // namespace Dune::BlockLab

#endif
