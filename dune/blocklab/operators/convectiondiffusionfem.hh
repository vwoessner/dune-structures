#ifndef DUNE_BLOCKLAB_OPERATORS_CONVECTIONDIFFUSIONFEM_HH
#define DUNE_BLOCKLAB_OPERATORS_CONVECTIONDIFFUSIONFEM_HH

/** Define a solver block that adds PDELab's convection diffusion
 *  FEM operator to the block solver parameter system.
 */

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/construction/muparser.hh>
#include<dune/blocklab/operators/virtualinterface.hh>
#include<dune/common/parametertree.hh>
#include<dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

#include<memory>
#include<string>


namespace Dune::BlockLab {


  /** A parameter class that fulfills the interface required by PDELab's
   *  convection diffusion operator, but parser the actual implementation
   *  from the inifile at runtime
   */
  template<typename GV, typename RF>
  class ConvectionDiffusionRuntimeParameters
  {
    public:
    using BCType = Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type;
    using Traits = Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF>;
    using RFT = typename Traits::RangeFieldType;
    using RT = typename Traits::RangeType;
    using E = typename Traits::ElementType;
    using ED = typename Traits::DomainType;
    using I = typename Traits::IntersectionType;
    using ID = typename Traits::IntersectionDomainType;
    using PT = typename Traits::PermTensorType;

    ConvectionDiffusionRuntimeParameters(const Dune::ParameterTree& config)
      : bctype_func(muparser_callable<BCType(I, ID)>(config.get<std::string>("bctype", "1")))
      , c_func(muparser_callable<RFT(E, ED)>(config.get<std::string>("c", "0.0")))
      , f_func(muparser_callable<RFT(E, ED)>(config.get<std::string>("f", "0.0")))
      , j_func(muparser_callable<RFT(I, ID)>(config.get<std::string>("j", "0.0")))
      , o_func(muparser_callable<RFT(I, ID)>(config.get<std::string>("o", "0.0")))
    {}

    // we cannot guarantee this statically, so we pay the price of always evaluating
    // at every quadrature point
    static constexpr bool permeabilityIsConstantPerCell()
    {
      return false;
    }

    // TODO Implement this!
    PT A (const E& e, const ED& x) const
    {
      PT I;
      for (std::size_t i=0; i<Traits::dimDomain; i++)
        for (std::size_t j=0; j<Traits::dimDomain; j++)
          I[i][j] = (i==j) ? 1 : 0;
      return I;
    }

    // TODO Implement this!
    RT b (const E& e, const ED& x) const
    {
      return RT(0.0);
    }

    //! reaction term
    RFT c (const E& e, const ED& x) const
    {
      return c_func(e, x);
    }

    //! source term
    RFT f (const E& e, const ED& x) const
    {
      return f_func(e, x);
    }

    //! boundary condition type function
    BCType bctype (const I& is, const ID& x) const
    {
      return bctype_func(is, x);
    }

    //! Neumann boundary condition
    RFT j (const I& is, const ID& x) const
    {
      return j_func(is, x);
    }

    //! outflow boundary condition
    RFT o (const I& is, const ID& x) const
    {
      return o_func(is, x);
    }

    private:
    std::function<BCType(I, ID)> bctype_func;
    std::function<RFT(E, ED)> c_func;
    std::function<RFT(E, ED)> f_func;
    std::function<RFT(I, ID)> j_func;
    std::function<RFT(I, ID)> o_func;
  };


  template<typename P, typename V, std::size_t i>
  class ConvectionDiffusionFEMBlock
    : public BlockBase<P, V, i>
  {
    public:
    using Traits = typename BlockBase<P, V, i>::Traits;
    using ParameterInterface = ConvectionDiffusionRuntimeParameters<typename Traits::EntitySet, typename Traits::Range>;
    using AbstractOperator = AbstractLocalOperatorInterface<typename Traits::GridFunctionSpace>;
    using LocalOperator = Dune::PDELab::ConvectionDiffusionFEM<ParameterInterface, typename Traits::FiniteElementMap>;
    using WrappedOperator = VirtualizedLocalOperator<LocalOperator, typename Traits::GridFunctionSpace>;

    template<typename Context>
    ConvectionDiffusionFEMBlock(Context&, const Dune::ParameterTree& config)
      : params(std::make_shared<ParameterInterface>(config))
    {}

    virtual ~ConvectionDiffusionFEMBlock() = default;

    virtual void setup() override
    {
      auto lop = std::make_shared<LocalOperator>(*params);
      auto wrapped = std::make_shared<WrappedOperator>(lop);
      this->solver->template introduce_parameter<std::shared_ptr<AbstractOperator>>("convectiondiffusionfem_operator", wrapped);
    }

    private:
    std::shared_ptr<ParameterInterface> params;
  };

} // namespace Dune::BlockLab

#endif
