FROM registry.dune-project.org/docker/ci/dune:git

ARG DUNECI_PARALLEL=1
ARG DUNECI_CMAKE_FLAGS="-DDUNE_PYTHON_VIRTUALENV_SETUP=True"

# Base requirements
USER root
RUN apt-get update \
    && apt-get -y upgrade \
    && apt-get -y install gmsh libyaml-cpp-dev python3-pip python3-venv \
    && apt-get clean
USER duneci

# Configure dune-venv
WORKDIR /duneci/modules
RUN cmake -DDUNE_PYTHON_VIRTUALENV_SETUP=True dune-common/build-cmake

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
