#ifndef DUNE_BLOCKLAB_CONSTRUCTION_CONTEXT_HH
#define DUNE_BLOCKLAB_CONSTRUCTION_CONTEXT_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/blocklab/solver.hh>
#include<dune/blocklab/utilities/enumerate.hh>
#include<dune/blocklab/utilities/stringsplit.hh>
#include<dune/common/hybridutilities.hh>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertree.hh>

#include<functional>
#include<memory>
#include<tuple>

namespace Dune::BlockLab {

  template<typename P, typename V>
  class ConstructionContext
  {
    public:
    ConstructionContext(Dune::MPIHelper& helper,
			const Dune::ParameterTree& rootconfig)
      : mpihelper(helper)
      , rootconfig(rootconfig)
    {}

    // The data members that are available for the construction of blocks
    Dune::MPIHelper& mpihelper;
    Dune::ParameterTree rootconfig;

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
      // Read the vector->index mapping from the configuration.
      // As soon as the interfaces for that are ready it should
      // be deduced from the vector providers!
      for (auto [i, v] : enumerate(string_split(config.get<std::string>("vectors"))))
	vector_names[v] = i;

      auto solver = std::make_shared<BlockSolver<P, V>>();
      solver->add(std::make_shared<ParentBlockBase<P, V>>(*this, config));
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
          auto vector = config.get<std::string>("vector");
          auto it = vector_names.find(vector);
          if (it == vector_names.end())
            DUNE_THROW(Dune::Exception, "The specified vector name does not exist");
          i = it->second;
        }

        identifier = "_" + identifier + "_" + std::to_string(i);
      }

      return mapping[identifier](*this, config);
    }

    private:
    // The construction function mapping
    std::map<std::string, std::function<std::shared_ptr<AbstractBlockBase<P, V>>(ConstructionContext<P, V>&, const Dune::ParameterTree&)>> mapping;
    std::map<std::string, std::size_t> vector_names;
  };

} // namespace Dune::BlockLab

#endif
