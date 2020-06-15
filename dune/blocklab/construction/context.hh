#ifndef DUNE_BLOCKLAB_CONSTRUCTION_CONTEXT_HH
#define DUNE_BLOCKLAB_CONSTRUCTION_CONTEXT_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/solver.hh>
#include<dune/blocklab/utilities/enumerate.hh>
#include<dune/blocklab/utilities/stringsplit.hh>
#include<dune/blocklab/utilities/tuplecat.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/hybridutilities.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertree.hh>

#include<array>
#include<functional>
#include<memory>
#include<tuple>

namespace Dune::BlockLab {

  template<typename UserParameters, typename... VectorProviders>
  class ConstructionContext
  {
    // Construct a tuple of all vector types in the given providers
    using V = std::tuple<typename VectorProviders::Vector...>;

    // and a tuple of all parameter types
    using P = tuple_cat_t<UserParameters, typename VectorProviders::Parameter...>;

    public:

    ConstructionContext(Dune::MPIHelper& helper,
			const Dune::ParameterTree& rootconfig,
			std::shared_ptr<VectorProviders>... providers)
      : mpihelper(helper)
      , rootconfig(rootconfig)
      , providers{providers...}
    {
      std::size_t i = 0;
      (this->vector_names.insert({{ providers->getName(), i++ }}) , ...);

      if (vector_names.size() != sizeof...(VectorProviders))
	DUNE_THROW(Dune::Exception, "Vector name was not unique throughout given vector providers");
    }

    // The following getter methods are accessible during block construction
    // through this Context object!
    const Dune::MPIHelper& getMPIHelper() const
    {
      return mpihelper;
    }

    const Dune::ParameterTree& getRootConfig() const
    {
      return rootconfig;
    }

    // Get a pointer to the solver object that is about to be constructed.
    // I am more than unhappy with this being in the interface because the
    // solver interface cannot be actually used during construction. But it
    // happens that some parts (Muparser callbacks) need the pointer to the
    // solver during construction! I am now handing out a const version,
    // which should offer some protection against misuse
    std::shared_ptr<const BlockSolver<P, V>> getSolver()
    {
      return solver;
    }

    template<typename Func>
    void registerBlock(std::string identifier, Func&& func)
    {
      mapping[identifier] = std::forward<Func>(func);
    }

    template<template<class, class, std::size_t> typename Step>
    void registerBlock(std::string identifier)
    {
      constexpr std::size_t vectors = Step<P, V, 0>::vectors;

      Dune::Hybrid::forEach(Dune::Hybrid::integralRange(std::integral_constant<std::size_t, 0>{},
                                                        std::integral_constant<std::size_t, std::tuple_size<V>::value - vectors + 1>{}),
			    [this, identifier](auto i)
			    {
                              auto subident = "_" + identifier + "_" + std::to_string(i);
                              this->registerBlock(subident,
						  [i](auto& ctx, const auto& p)
						  {
                        	                    return std::make_shared<Step<P, V, i>>(ctx, p);
						  });
			    });
    }

    template<template<class, class> typename Step>
    void registerBlock(std::string identifier)
    {
      registerBlock(identifier,
		    [](auto& ctx, const auto& p)
		    {
	              return std::make_shared<Step<P, V>>(ctx, p);
		    });
    }

    std::shared_ptr<BlockSolver<P, V>> constructSolver(const Dune::ParameterTree& config)
    {
      if(!solver)
      {
	solver = std::make_shared<BlockSolver<P, V>>(
	  std::apply([](auto... p){ return std::make_tuple(p->getVector()...); }, providers),
	  std::apply([](auto... p){ return std::make_tuple(p->getConstraintsContainer()...); }, providers)
	);

	solver->add(std::make_shared<ParentBlockBase<P, V>>(*this, config));
      }
      return solver;
    }

    std::shared_ptr<AbstractBlockBase<P, V>> constructBlock(std::string stepname, const Dune::ParameterTree& config)
    {
      auto identifier = config.get<std::string>("type", stepname);

      // If this identifier does not exist, we are treating it as a vector identifier
      if (mapping.find(identifier) == mapping.end())
      {
        // Find the correct vector index
	std::size_t i = 0;
	if (config.hasKey("vector"))
        {
          auto vector = config.get<std::string>("vector", "solution");
          auto it = vector_names.find(vector);
          if (it == vector_names.end())
            DUNE_THROW(Dune::Exception, "The specified vector name does not exist");
          i = it->second;
        }

        identifier = "_" + identifier + "_" + std::to_string(i);
      }

      std::shared_ptr<AbstractBlockBase<P, V>> block;
      try {
        block = mapping[identifier](*this, config);
      }
      catch(std::bad_function_call&)
      {
        DUNE_THROW(Dune::Exception, "Looking for a block of type '" + identifier + "' failed!");
      }

      return block;
    }

    private:
    // The construction function mapping
    std::map<std::string, std::function<std::shared_ptr<AbstractBlockBase<P, V>>(ConstructionContext<UserParameters, VectorProviders...>&, const Dune::ParameterTree&)>> mapping;
    std::map<std::string, std::size_t> vector_names;

    // The vector providers - maybe we can achieve that all information
    // is extracted already and not store this!
    std::tuple<std::shared_ptr<VectorProviders>...> providers;

    // The objects that are available during grid construction through above getter
    Dune::MPIHelper& mpihelper;
    Dune::ParameterTree rootconfig;

    // The solver class. Actually this should not be exposed during construction
    // because only a very limited subset of its interface is valid
    std::shared_ptr<BlockSolver<P, V>> solver;
  };

} // namespace Dune::BlockLab

#endif
