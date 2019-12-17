add_executable(cubeapp cubeapp.cc)
add_dependencies(cubeapp operators)
dune_symlink_to_source_files(FILES iterative_compress.ini linearized_compress.ini)
