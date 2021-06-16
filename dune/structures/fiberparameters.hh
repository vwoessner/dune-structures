#ifndef DUNE_STRUCTURES_FIBERPARAMETERS_HH
#define DUNE_STRUCTURES_FIBERPARAMETERS_HH

#include <yaml-cpp/yaml.h>

/// A storage for material parameters of a fiber
template<typename Range = double>
struct FiberParameters
{
  Range radius;
  Range youngs_modulus;
  Range prestress;

  /// Default constructor
  FiberParameters() = default;
  /// Default destructor
  ~FiberParameters() = default;

  /// Constructor from explicit input values
  FiberParameters(const Range radius_,
                  const Range youngs_modulus_,
                  const Range prestress_)
    : radius(radius_)
    , youngs_modulus(youngs_modulus_)
    , prestress(prestress_)
  {
  }

  /// Constructor from a YAML config node
  FiberParameters(const YAML::Node& config)
    : radius(config["radius"].as<Range>())
    , youngs_modulus(config["youngs_modulus"].as<Range>())
    , prestress(config["prestress"].as<Range>())
  {
  }
};

#endif // DUNE_STRUCTURES_FIBERPARAMETERS_HH
