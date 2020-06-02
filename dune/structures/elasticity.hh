#ifndef DUNE_STRUCTURES_ELASTICITY_HH
#define DUNE_STRUCTURES_ELASTICITY_HH

#include<dune/pdelab.hh>

#include<iostream>
#include<memory>
#include<tuple>

// The 3D operators
#include"operators/elasticity_operator.hh"
#include"operators/quasistatic_mass_operator.hh"
#include"operators/elasticity_p2_operator.hh"
#include"operators/quasistatic_mass_p2_operator.hh"

// The 2D operators
#include"operators/elasticity_2d_operator.hh"
#include"operators/quasistatic_mass_2d_operator.hh"
#include"operators/elasticity_2d_p2_operator.hh"
#include"operators/quasistatic_mass_2d_p2_operator.hh"


template<typename GFS, int dim, typename FEM>
struct OperatorSwitchImpl
{};

template<typename GFS, int dim>
using OperatorSwitch = OperatorSwitchImpl<GFS, dim, typename GFS::template Child<0>::Type::Traits::FiniteElementMap>;

template<typename GFS, typename ES, typename R>
struct OperatorSwitchImpl<GFS, 3, Dune::PDELab::PkLocalFiniteElementMap<ES, double, R, 1>>
{
  using Elasticity = ElasticityOperator<GFS, GFS, GFS, GFS>;
  using Mass = QuasiStaticMassOperator<GFS, GFS>;
};

template<typename GFS, typename ES, typename R>
struct OperatorSwitchImpl<GFS, 3, Dune::PDELab::PkLocalFiniteElementMap<ES, double, R, 2>>
{
  using Elasticity = ElasticityP2Operator<GFS, GFS, GFS, GFS>;
  using Mass = QuasiStaticMassP2Operator<GFS, GFS>;
};

template<typename GFS, typename ES, typename R>
struct OperatorSwitchImpl<GFS, 2, Dune::PDELab::PkLocalFiniteElementMap<ES, double, R, 1>>
{
  using Elasticity = Elasticity2DOperator<GFS, GFS, GFS, GFS>;
  using Mass = QuasiStaticMass2DOperator<GFS, GFS>;
};

template<typename GFS, typename ES, typename R>
struct OperatorSwitchImpl<GFS, 2, Dune::PDELab::PkLocalFiniteElementMap<ES, double, R, 2>>
{
  using Elasticity = Elasticity2DP2Operator<GFS, GFS, GFS, GFS>;
  using Mass = QuasiStaticMass2DP2Operator<GFS, GFS>;
};


template<int degree, typename ES, typename RangeType = double>
auto elasticity_setup(ES es)
{
  constexpr int dim = ES::dimension;

  // Set up finite element maps...
  using FEM = Dune::PDELab::PkLocalFiniteElementMap<ES, double, RangeType, degree>;
  auto fem = std::make_shared<FEM>(es);

  // Set up grid function spaces...
  using VB = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
  using CASS = Dune::PDELab::ConformingDirichletConstraints;
  using GFS = Dune::PDELab::VectorGridFunctionSpace<ES, FEM, dim, VB, VB, CASS>;
  auto gfs = std::make_shared<GFS>(es, fem);
  gfs->name("displacement");
  gfs->update();
  std::cout << "Set up a grid function space with " << gfs->size() << " dofs!" << std::endl;

  // Setting up constraints container
  using CC = typename GFS::template ConstraintsContainer<RangeType>::Type;
  auto cc = std::make_shared<CC>();
  cc->clear();

  // Setting up container
  using V = Dune::PDELab::Backend::Vector<GFS, double>;
  auto x = std::make_shared<V>(gfs);

  // And two more containers for body force and traction fields. One might consider
  // interpolation with lower order here.
  auto force = std::make_shared<V>(gfs, 0.0);
  auto traction = std::make_shared<V>(gfs, 0.0);

  return std::make_tuple(x, cc, force, traction);
}


template<int degree, typename ES, typename RangeType = double>
auto elastodynamics_setup(ES es)
{
  // Set up finite element maps...
  using FEM = Dune::PDELab::PkLocalFiniteElementMap<ES, double, RangeType, degree>;
  auto fem = std::make_shared<FEM>(es);

  // Set up grid function spaces...
  using VB = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
  using CASS = Dune::PDELab::ConformingDirichletConstraints;
  using GFS = Dune::PDELab::VectorGridFunctionSpace<ES, FEM, 3, VB, VB, CASS>;
  using PGFS = Dune::PDELab::CompositeGridFunctionSpace<VB, Dune::PDELab::LexicographicOrderingTag, GFS, GFS>;

  auto gfs1 = std::make_shared<GFS>(es, fem);
  auto gfs2 = std::make_shared<GFS>(es, fem);
  auto pgfs = std::make_shared<PGFS>(gfs1, gfs2);
  gfs1->name("displacement");
  gfs2->name("auxiliary");
  pgfs->update();
  std::cout << "Set up a grid function space with " << pgfs->size() << " dofs!" << std::endl;

  // Setting up constraints container
  using CC = typename PGFS::template ConstraintsContainer<RangeType>::Type;
  auto cc = std::make_shared<CC>();
  cc->clear();

  // Setting up container
  using V = Dune::PDELab::Backend::Vector<PGFS, double>;
  auto x = std::make_shared<V>(pgfs);

  return std::make_pair(x, cc);
}


#endif
