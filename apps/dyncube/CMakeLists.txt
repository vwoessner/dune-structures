add_executable(dyncubeapp dyncubeapp.cc)
add_dependencies(dyncubeapp operators)
dune_symlink_to_source_files(FILES dyncubeapp.ini)
