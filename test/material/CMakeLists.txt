dune_add_test(SOURCES testmaterial.cc)

dune_structures_mesh_generation(TARGET testmaterial
                                CONFIG materialgrid.yml)
