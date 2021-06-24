dune_add_test(SOURCES test-parameter-assign.cc)
dune_add_test(SOURCES test-fiber-intersection.cc)

configure_file(cube.msh cube.msh COPYONLY)
add_test(NAME run_gmsh
         COMMAND gmsh2d_p2 run ${CMAKE_CURRENT_SOURCE_DIR}/cube.yml
)
