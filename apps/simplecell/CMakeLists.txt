add_executable(simplecell simplecell.cc)
add_dependencies(simplecell operators)

dune_structures_mesh_generation(TARGET simplecell
                                CONFIG simplecell.ini
                                )

dune_symlink_to_source_files(FILES simplecell.ini)
