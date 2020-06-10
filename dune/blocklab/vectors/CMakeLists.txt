add_subdirectory(test)

install(FILES concept.hh
              pkfem.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/blocklab/vectors)
