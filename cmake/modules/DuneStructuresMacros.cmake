# Make the code generator sources of dune-structures known
# to dune-codegen to correctly retrigger code generation upon changes
# in the python code.

file(GLOB_RECURSE _ADDITIONAL_SOURCES "${CMAKE_SOURCE_DIR}/python/*.py")
set(DUNE_CODEGEN_ADDITIONAL_PYTHON_SOURCES ${DUNE_CODEGEN_ADDITIONAL_PYTHON_SOURCES} ${_ADDITIONAL_SOURCES})


function(dune_structures_mesh_generation)
  include(CMakeParseArguments)
  cmake_parse_arguments(MESHGEN "" "TARGET;CONFIG" "" ${ARGN})

  # Add a build rule for the mesh file
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${MESHGEN_TARGET}.msh
                     COMMAND ${CMAKE_BINARY_DIR}/run-in-dune-env generate_cell_mesh ${CMAKE_CURRENT_SOURCE_DIR}/${MESHGEN_CONFIG} ${CMAKE_CURRENT_BINARY_DIR}/${MESHGEN_TARGET}.msh
                     DEPENDS ${CMAKE_SOURCE_DIR}/python/dune/structures/gmsh.py ${CMAKE_CURRENT_SOURCE_DIR}/${MESHGEN_CONFIG}
                     )

  # Have the executable depend on the meshfile
  target_sources(${MESHGEN_TARGET} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/${MESHGEN_TARGET}.msh)
endfunction()