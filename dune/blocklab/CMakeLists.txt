add_subdirectory(blocks)
add_subdirectory(construction)
add_subdirectory(grids)
add_subdirectory(operators)
add_subdirectory(test)
add_subdirectory(utilities)
add_subdirectory(vectors)

install(FILES solver.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/blocklab
        )
