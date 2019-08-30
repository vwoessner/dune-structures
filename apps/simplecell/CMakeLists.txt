dune_add_generated_executable(SOURCE simplecell.cc
                              INIFILE simplecell.ini
                              UFLFILE simplecell.ufl
                              TARGET simplecell)

dune_structures_mesh_generation(TARGET simplecell
                                CONFIG simplecell.ini
                                )

dune_symlink_to_source_files(FILES simplecell.ini)
