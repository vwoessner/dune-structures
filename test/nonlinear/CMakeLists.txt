dune_add_formcompiler_system_test(SOURCE testnonlinear.cc
                                  INIFILE nonlinear.ini
                                  UFLFILE nonlinear.ufl
                                  BASENAME testnonlinear)

dune_symlink_to_source_files(FILES nonlinear.ini)