FROM registry.dune-project.org/docker/ci/dune:git

ARG DUNECI_PARALLEL=1
ARG DUNECI_CMAKE_FLAGS="-DDUNE_PYTHON_VIRTUALENV_SETUP=True"

# Base requirements
USER root
RUN apt-get update \
    && apt-get -y upgrade \
    && apt-get -y install gmsh libgl1-mesa-dev libyaml-cpp-dev python3-pip python3-venv \
    && apt-get clean
USER duneci

# Configure dune-venv
WORKDIR /duneci/modules
RUN cmake -DDUNE_PYTHON_VIRTUALENV_SETUP=True dune-common/build-cmake

# BIG: Build VTK wheels
WORKDIR /duneci/vtk
RUN git clone https://github.com/Kitware/VTK -b release \
    && mkdir VTK/build && cd VTK/build \
    && cmake -DCMAKE_BUILD_TYPE=Release -DVTK_BUILD_TESTING=OFF \
        -DVTK_BUILD_DOCUMENTATION=OFF -DVTK_BUILD_EXAMPLES=OFF -DVTK_DATA_EXCLUDE_FROM_ALL:BOOL=ON \
        -DVTK_MODULE_ENABLE_VTK_PythonInterpreter:STRING=NO \
        -DVTK_WHEEL_BUILD=ON -DVTK_PYTHON_VERSION=3 -DVTK_WRAP_PYTHON=ON \
        ../ \
    && make -j$DUNECI_PARALLEL \
    && /duneci/modules/dune-common/build-cmake/dune-env/bin/python3 -m pip install setuptools wheel \
    && /duneci/modules/dune-common/build-cmake/dune-env/bin/python3 setup.py bdist_wheel \
    && /duneci/modules/dune-common/build-cmake/dune-env/bin/python3 -m pip install dist/vtk-*.whl \
    && cd /duneci && rm -rf vtk/
WORKDIR /duneci/modules

# alugrid, testtools, pdelab
RUN duneci-install-module --recursive https://gitlab.dune-project.org/extensions/dune-alugrid.git
RUN duneci-install-module --recursive https://gitlab.dune-project.org/quality/dune-testtools.git
RUN duneci-install-module --recursive https://gitlab.dune-project.org/pdelab/dune-pdelab.git

# Build dune-codegen after applying patches
WORKDIR /duneci/modules
RUN git clone --recurse-submodules https://gitlab.dune-project.org/extensions/dune-codegen.git
RUN cd /duneci/modules/dune-codegen && bash patches/apply_patches.sh && cd /duneci/modules/ \
    && dunecontrol --only=dune-codegen --opts="/duneci/dune.opts" all

# Build dune-blocklab after adding frontend requirements
WORKDIR /duneci/modules
RUN git clone --recurse-submodules https://gitlab.dune-project.org/dominic/dune-blocklab.git
RUN python3 -m pip install -r dune-blocklab/frontend/requirements.txt
RUN dunecontrol --only=dune-blocklab --opts="/duneci/dune.opts" all
