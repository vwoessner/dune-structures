ARG TOOLCHAIN="gcc-9-17"
FROM registry.dune-project.org/dominic/dune-blocklab/frontend:${TOOLCHAIN}

RUN duneci-install-module --recursive https://gitlab.dune-project.org/dominic/dune-structures.git
