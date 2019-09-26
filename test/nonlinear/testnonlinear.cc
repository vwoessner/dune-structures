#include"config.h"

#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/uggrid.hh>
#include<dune/pdelab.hh>
#include<dune/structures/material.hh>
#include<dune/structures/onetoone.hh>
#include<dune/structures/vonmises.hh>
#include<dune/structures/visualization.hh>
#include<dune/testtools/gridconstruction.hh>

#include<vector>

#include"elasticity_operator.hh"

int main(int argc, char** argv)
{
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  // Parse the ini file
  Dune::ParameterTree params;
  Dune::ParameterTreeParser::readINITree("nonlinear.ini", params);

  using RangeType = double;

  // Construct the grid and loadbalance it
  using GridType = Dune::UGGrid<3>;
  using GV = GridType::LeafGridView;
  using ES = Dune::PDELab::NonOverlappingEntitySet<GV>;
  using DF = GridType::ctype;

  Dune::FieldVector<double, 3> ll(0.0);
  Dune::FieldVector<double, 3> ur(1.0);
  auto N = params.get<std::array<unsigned int, 3>>("grid.N", {10, 10, 10});
  std::cout << "N: " << N[0] << " " << N[1] << " " << N[2] << std::endl;

  auto grid = Dune::StructuredGridFactory<GridType>::createSimplexGrid(ll, ur, N);
  auto physical = std::make_shared<std::vector<int>>(grid->size(0), 0);
  GV gv = grid->leafGridView();
  ES es(grid->leafGridView());

  // Set up finite element maps...
  using FEM = Dune::PDELab::PkLocalFiniteElementMap<ES, DF, RangeType, 1>;
  FEM fem(es);

  // Set up grid function spaces...
  using VB = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
  using CASS = Dune::PDELab::ConformingDirichletConstraints;
  using GFS = Dune::PDELab::VectorGridFunctionSpace<ES, FEM, 3, VB, VB, CASS>;
  GFS gfs(es, fem);
  gfs.name("displacement");
  gfs.update();
  std::cout << "Set up a grid function space with " << gfs.size() << " dofs!" << std::endl;

  // Setting up constraints container
  using CC = GFS::ConstraintsContainer<RangeType>::Type;
  CC cc;
  cc.clear();
  auto bctype = [&](const auto& x){ return (x[2] < 1e-08) || (x[2] > 1.0 - 1e-8); };
  auto bctype_f = Dune::PDELab::makeBoundaryConditionFromCallable(es, bctype);
  Dune::PDELab::CompositeConstraintsParameters comp_bctype(bctype_f, bctype_f, bctype_f);
  Dune::PDELab::constraints(comp_bctype, gfs, cc);
  std::cout << "Set up a constraints container with " << cc.size() << " dofs!" << std::endl;

  // Instantiate the material class
  auto material = parse_material<RangeType>(es, physical, params.sub("material"));

  // Setting up grid operator
  using LOP = ElasticityOperator<GFS, GFS>;
  using MB = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  using GO = Dune::PDELab::GridOperator<GFS, GFS, LOP, MB, DF, RangeType, RangeType, CC, CC>;
  LOP lop(gfs, gfs, params, material);
  MB mb(21);
  GO go(gfs, cc, gfs, cc, lop, mb);

  // Setting up container
  using V = Dune::PDELab::Backend::Vector<GFS, DF>;
  V x(gfs);
  auto compression = params.get<double>("model.compression");
  auto dirichlet = [&](const auto& x){ return (compression - 1.0) * x[2]; };
  auto zero = [](const auto& x){ return 0.0; };
  auto zero_f = Dune::PDELab::makeGridFunctionFromCallable(es, zero);
  auto dirichlet_f = Dune::PDELab::makeGridFunctionFromCallable(es, dirichlet);
  Dune::PDELab::CompositeGridFunction comp_dirichlet(zero_f, zero_f, dirichlet_f);
  Dune::PDELab::interpolate(comp_dirichlet, gfs, x);

  // Assert sanity of the initial condition
  Dune::PDELab::VectorDiscreteGridFunction<GFS, V> gf(gfs, x);
  if (!is_onetoone(gf))
    DUNE_THROW(Dune::Exception, "Self-intersecting displacement field detected");

  // Set up the solver...
  using LS = Dune::PDELab::ISTLBackend_SEQ_UMFPack;
//  using LS = Dune::PDELab::ISTLBackend_NOVLP_BCGS_AMG_SSOR<GO>;
//  using LS = Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<GO>;

  using NLP = Dune::PDELab::Newton<GO, LS, V>;
//  LS ls(go, 5000, 3);
  LS ls;
  NLP nlp(go, x, ls);
  nlp.setParameters(params.sub("newton"));

  try {
    nlp.apply();
  }
  catch (Dune::PDELab::NewtonError& e)
  {
    std::cout << "Encountered a fatal Newton error. Diagnosing..." << std::endl;
    if(!is_onetoone(gf))
    {
      is_onetoone(gf, true);
      std::cout << "Self-intersecting displacement field detected!" << std::endl;
    }
    else
      std::cout << "Self-intersecting was not the problem!" << std::endl;
    return 1;
  }

  // A grid function for the stress
  VonMisesStressGridFunction stress(x, material);

  // Interpolate the stress into a grid function
  using P0FEM = Dune::PDELab::P0LocalFiniteElementMap<DF, RangeType, 3>;
  P0FEM p0fem(Dune::GeometryTypes::simplex(3));
  using P0GFS = Dune::PDELab::GridFunctionSpace<ES, P0FEM, CASS, VB>;
  P0GFS p0gfs(es, p0fem);
  p0gfs.name("vonmises");
  using SV = Dune::PDELab::Backend::Vector<P0GFS, DF>;
  SV stress_container(p0gfs);
  Dune::PDELab::interpolate(stress, p0gfs, stress_container);

  // Visualize the stress grid function
  Dune::VTKWriter vtkwriter(es.gridView());
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, x);
  Dune::PDELab::addSolutionToVTKWriter(vtkwriter, p0gfs, stress_container);
  vtkwriter.addCellData(*physical, "gmshPhysical");
  write_rankdata(vtkwriter, helper, es.gridView());
  vtkwriter.write("output", Dune::VTK::ascii);

  return 0;
}
