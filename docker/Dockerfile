FROM debian:buster
MAINTAINER Dominic Kempf (dominic.kempf@iwr.uni-heidelberg.de)

RUN adduser --disabled-password --home /structure --uid 50000 structure

# Debian packages
RUN export DEBIAN_FRONTEND=noninteractive; \
  apt-get update && apt-get dist-upgrade --no-install-recommends --yes \
  && apt-get install --no-install-recommends --yes \
    build-essential \
    ca-certificates \
    cmake \
    gfortran \
    git \
    gmsh \
    libisl-dev \
    libmuparser-dev \
    libopenmpi-dev \
    libsuitesparse-dev \
    openmpi-bin \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

# Clone all the necessary Dune modules
RUN cd /structure \
 && git clone https://gitlab.dune-project.org/core/dune-common.git \
 && git clone https://gitlab.dune-project.org/core/dune-grid.git \
 && git clone https://gitlab.dune-project.org/core/dune-geometry.git \
 && git clone https://gitlab.dune-project.org/core/dune-istl.git \
 && git clone https://gitlab.dune-project.org/core/dune-localfunctions.git \
 && git clone https://gitlab.dune-project.org/staging/dune-typetree.git \
 && git clone https://gitlab.dune-project.org/staging/dune-functions.git \
 && git clone https://gitlab.dune-project.org/staging/dune-uggrid.git \
 && git clone https://gitlab.dune-project.org/pdelab/dune-pdelab.git \
 && git clone https://gitlab.dune-project.org/quality/dune-testtools.git \
 && git clone https://gitlab.dune-project.org/extensions/dune-alugrid.git \
 && git clone --recursive https://gitlab.dune-project.org/extensions/dune-codegen.git

COPY config.opts /structure/

RUN /structure/dune-common/bin/dunecontrol --opts=/structure/config.opts --module=dune-codegen --builddir=/structure/build all

RUN git clone --recursive https://gitlab.dune-project.org/dominic/dune-structures.git
