add_subdirectory(test)

install(FILES convectiondiffusionfem.hh
              virtualinterface.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/blocklab/operators)