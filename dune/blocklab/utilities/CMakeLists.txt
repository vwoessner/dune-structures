add_subdirectory(test)

install(FILES enumerate.hh
              stringsplit.hh
              tuplecat.hh
              uniquevariant.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/blocklab/utilities)