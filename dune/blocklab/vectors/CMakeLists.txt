add_subdirectory(test)

install(FILES composite.hh
              concept.hh
              pkfem.hh
              power.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/blocklab/vectors)
