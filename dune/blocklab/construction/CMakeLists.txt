add_subdirectory(test)

install(FILES callabletree.hh
              context.hh
              muparser.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/blocklab/construction)