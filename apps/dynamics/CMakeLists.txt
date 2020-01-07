add_executable(dynamicsapp dynamicsapp.cc)
add_dependencies(dynamicsapp operators)

dune_symlink_to_source_files(FILES linear_dyncube.ini
                                   neohookean_dyncube.ini
                                   stvenant_dyncube.ini
                             )
