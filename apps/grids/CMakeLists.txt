add_executable(testgrid testgrid.cc)

dune_structures_mesh_generation(TARGET testgrid CONFIG testgrid.yml)
