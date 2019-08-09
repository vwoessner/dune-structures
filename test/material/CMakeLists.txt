dune_add_test(SOURCES testmaterial.cc)

dune_structures_mesh_generation(TARGET testmaterial
                                CONFIG materialgrid.yml)

dune_symlink_to_source_files(FILES material.ini)
