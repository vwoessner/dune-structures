add_subdirectory(test)

install(FILES adaptivity.hh
              blockbase.hh
              blocktraits.hh
              constraints.hh
              control.hh
              filelogger.hh
              instationary.hh
              linearsolver.hh
              newton.hh
              parameter.hh
              probe.hh
              variation.hh
              visualization.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/blocklab/blocks)