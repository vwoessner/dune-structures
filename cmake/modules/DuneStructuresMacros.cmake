# Make the code generator sources of dune-structures known
# to dune-codegen to correctly retrigger code generation upon changes
# in the python code.

file(GLOB_RECURSE _ADDITIONAL_SOURCES "${CMAKE_SOURCE_DIR}/python/dune/structures/*.py")
set(DUNE_CODEGEN_ADDITIONAL_PYTHON_SOURCES ${DUNE_CODEGEN_ADDITIONAL_PYTHON_SOURCES} ${_ADDITIONAL_SOURCES})

# Search for additional packages needed by dune-structures
find_package(ParMETIS 4)
find_package(Gmsh)

include(pandocology)
