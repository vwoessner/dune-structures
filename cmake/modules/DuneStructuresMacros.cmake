# Make the code generator sources of dune-structures known
# to dune-codegen to correctly retrigger code generation upon changes
# in the python code.

file(GLOB_RECURSE _ADDITIONAL_SOURCES "${CMAKE_SOURCE_DIR}/python/*.py")
set(DUNE_CODEGEN_ADDITIONAL_PYTHON_SOURCES ${DUNE_CODEGEN_ADDITIONAL_PYTHON_SOURCES} ${_ADDITIONAL_SOURCES})
