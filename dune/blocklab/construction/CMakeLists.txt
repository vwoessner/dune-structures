add_subdirectory(test)

install(FILES callabletree.hh
              context.hh
              dataparser.hh
              muparser.hh
              registry.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/blocklab/construction)