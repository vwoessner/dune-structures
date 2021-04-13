# Make the code generator sources of dune-structures known
# to dune-codegen to correctly retrigger code generation upon changes
# in the python code.
set(DUNE_CODEGEN_ADDITIONAL_PYTHON_SOURCES
    ${DUNE_CODEGEN_ADDITIONAL_PYTHON_SOURCES}
    ${PROJECT_SOURCE_DIR}/python/dune/structures/__init__.py
    ${PROJECT_SOURCE_DIR}/python/dune/structures/adhesion.py
    ${PROJECT_SOURCE_DIR}/python/dune/structures/codegen.py
    ${PROJECT_SOURCE_DIR}/python/dune/structures/diffgeo.py
    ${PROJECT_SOURCE_DIR}/python/dune/structures/elasticity.py
    ${PROJECT_SOURCE_DIR}/python/dune/structures/gmsh.py
    ${PROJECT_SOURCE_DIR}/python/dune/structures/material.py
)

# Generate the operators with a dummy executable (we only need the headers)
dune_add_generated_executable(INIFILE operators.ini
                              UFLFILE operators.ufl
                              TARGET operators
                              SOURCE dummy.cc)
