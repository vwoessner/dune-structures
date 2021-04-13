set(APPS
    structured2d_p2
    structured3d_p1
    cell3d_p1
    cell3d_p2
)

# Add apps as executables, register with frontend, make depend on code generator
foreach(exe ${APPS})
    add_executable(${exe} ${exe}.cc)
    blocklab_add_app_to_frontend(TARGET ${exe})
    add_dependencies(${exe} generation)
endforeach()

# Register tests
foreach(num 23_1 23_2 23_3 24_1 24_2)
    add_test(NAME run_hansbo_${num}
        COMMAND structured2d_p2 run ${PROJECT_SOURCE_DIR}/experiments/hansbo_reproduction/fig${num}.yml
)
endforeach()
