dune_add_generated_executable(SOURCE cubeapp.cc
                              INIFILE cubeapp.ini
                              UFLFILE cubeapp.ufl
                              TARGET cubeapp)

dune_symlink_to_source_files(FILES cubeapp.ini)
