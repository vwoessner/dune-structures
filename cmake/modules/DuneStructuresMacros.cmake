# Make the code generator sources of dune-structures known
# to dune-codegen to correctly retrigger code generation upon changes
# in the python code.

file(GLOB_RECURSE _ADDITIONAL_SOURCES "${CMAKE_SOURCE_DIR}/python/dune/structures/*.py")
set(DUNE_CODEGEN_ADDITIONAL_PYTHON_SOURCES ${DUNE_CODEGEN_ADDITIONAL_PYTHON_SOURCES} ${_ADDITIONAL_SOURCES})


# Search for additional packages needed by dune-structures
find_package(ParMETIS 4)
find_package(Gmsh)

# Implement automatic, CMake-triggered grid generation
function(dune_structures_mesh_generation)
  include(CMakeParseArguments)
  cmake_parse_arguments(MESHGEN "" "TARGET" "CONFIG" ${ARGN})

  # Gather pygmsh sources to retrigger mesh generation upon their changes
  file(GLOB_RECURSE PYGMSH_SOURCES "${CMAKE_SOURCE_DIR}/python/pygmsh/*.py")

  # Add a build rule for the mesh file
  foreach(conf ${MESHGEN_CONFIG})
    get_filename_component(cutname "${conf}" NAME_WE)
    add_custom_command(OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/${cutname}.msh"
                       COMMAND "${CMAKE_BINARY_DIR}/run-in-dune-env" generate_cell_mesh
                               "${CMAKE_CURRENT_SOURCE_DIR}/${conf}"
                               "${CMAKE_CURRENT_BINARY_DIR}/${cutname}.msh"
                               "${GMSH_EXECUTABLE}"
                       DEPENDS "${CMAKE_SOURCE_DIR}/python/dune/structures/gmsh.py"
                               "${CMAKE_CURRENT_SOURCE_DIR}/${conf}"
                               "${PYGMSH_SOURCES}"
                       )

    # Have the executable depend on the meshfile
    target_sources(${MESHGEN_TARGET} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/${cutname}.msh)
  endforeach()
endfunction()
