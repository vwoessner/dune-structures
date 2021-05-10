cmake_minimum_required(VERSION 3.16)
project(dune-structures CXX)

# Set up dune-common
if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

# Find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})

set(DUNE_PYTHON_INSTALL_EDITABLE 1)

#include the dune macros
set(DUNE_REENABLE_ADD_TEST TRUE)
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

dune_enable_all_packages()

add_subdirectory(doc)
add_subdirectory(operators)
add_subdirectory(apps)
add_subdirectory(cmake/modules)
add_subdirectory(dune)
add_subdirectory(python)
add_subdirectory(experiments)
add_subdirectory(test)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
