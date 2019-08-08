dune_add_test(SOURCES testgrid.cc)

dune_structures_mesh_generation(TARGET testgrid
                                CONFIG testgrid.yml)
