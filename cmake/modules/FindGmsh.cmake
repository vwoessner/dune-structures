# Check the availability of the Gmsh command.
#

find_program(GMSH_EXECUTABLE gmsh)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "Gmsh"
  DEFAULT_MSG
  GMSH_EXECUTABLE
)