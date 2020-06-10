add_subdirectory(test)

install(FILES capabilities.hh
              concept.hh
              gmsh.hh
              structured.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/blocklab/grids)