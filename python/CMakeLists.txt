# Make sure that this repository was cloned recursively
include(DuneSubmoduleChecker)
dune_check_submodule(PATH pygmsh INDICATOR setup.py)

# Install pygmsh manually in case we want to fiddle around with it
dune_python_install_package(PATH pygmsh)

# Install the Python package dune.structures into the Dune virtual environment
dune_python_install_package(PATH .)

# Add Python tests
dune_python_add_test(NAME python_tests
    INTERPRETER ${DUNE_PYTHON_VIRTUALENV_EXECUTABLE}
    MODULE pytest -v
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/dune/structures/test
)
