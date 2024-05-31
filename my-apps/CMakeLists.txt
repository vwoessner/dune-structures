set(APPS
    my_structured2d_p2
    my_cell2d_p2
    my_gmsh2d_p2
)

# Add apps as executables, register with frontend, make depend on code generator
foreach(exe ${APPS})
    add_executable(${exe} ${exe}.cc)
    add_dependencies(${exe} generation)
endforeach()

dune_symlink_to_source_files(FILES scheiwe.msh)
