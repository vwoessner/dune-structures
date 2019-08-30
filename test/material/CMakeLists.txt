dune_add_formcompiler_system_test(SOURCE testmaterial.cc
                                  INIFILE material.ini
                                  UFLFILE dummy.ufl
                                  BASENAME testmaterial)

dune_structures_mesh_generation(TARGET testmaterial
                                CONFIG material.ini)
