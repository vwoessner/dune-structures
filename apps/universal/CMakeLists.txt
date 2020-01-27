add_executable(universalapp universalapp.cc)
add_dependencies(universalapp operators)

dune_symlink_to_source_files(FILES fibrecube.ini
                                   iterative_cube_compress.ini
                                   linearized_cube_compress.ini
                                   quasistatic.ini
                                   robert.ini
                                   simplecell.ini
                                   twod_beam.ini
                             )
