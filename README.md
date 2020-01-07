This module contains the source code for the simulation of single cell mechanics and stress fibre formation executed in the STRUCTURES cluster of excellence at Heidelberg university. For questions write to [dominic.kempf@iwr.uni-heidelberg.de](mailto:dominic.kempf@iwr.uni-heidelberg.de).

# 	Installing the software

There is two principal ways to use this software: As dockerized applications or by building it from source. The former is recommended if you intend to *run* simulations, the latter is necessary if you want to *develop* your own simulators.

## Dockerized setup

Make sure that you have Docker installed on your machine.
You can do so e.g. by following [these instructions from Docker](https://docs.docker.com/install/linux/docker-ce/debian/#install-using-the-repository). Having done so, you can run the container interactively by doing:

```
docker run -ti registry.dune-project.org/dominic/dune-structures
```

More information that should be provided here:
* Entrypoints of the Docker container
* Mounting volumes into the container to retrieve visualization data

## Full installation

This is not an installation guide, but merely documents the specific requirements of dune-structures.
It assumes that you know how to build a Dune software stack.

dune-structures is a Dune module and as such has some external dependencies:

* CMake >= 3.1
* C++17-compliant compiler e.g. gcc >= 8.1
* MPI (e.g. Debian packages `libopenmpi-dev` and `openmpi-bin`)

dune-structures requires the following Dune modules:

* All the [Dune core modules](https://gitlab.dune-project.org/core)
* [dune-typetree](https://gitlab.dune-project.org/staging/dune-typetree.git)
* [dune-functions](https://gitlab.dune-project.org/staging/dune-functions.git)
* [dune-pdelab](https://gitlab.dune-project.org/pdelab/dune-pdelab.git)
* [dune-testtools](https://gitlab.dune-project.org/quality/dune-testtools.git)
* [dune-uggrid](https://gitlab.dune-project.org/staging/dune-uggrid.git)
* [dune-codegen](https://gitlab.dune-project.org/extensions/dune-codegen.git)

Note that dune-codegen requires special care in setup. It needs to
* be cloned with the `--recursive` flag
* You need to run `./patches/apply_patches.sh` from the module root directory after cloning

Currently, the latest versions (master) of all Dune modules are necessary.
Furthermore, dune-structures requires the following external packages:

* A Python interpreter >= 3.6 (e.g. Debian packages `python3-dev`)
* pip (e.g. the Debian package `python3-pip`)
* MuParser (e.g. the Debian package `libmuparser-dev`)
* Gmsh (e.g. Debian packages `gmsh`)

You should build your Dune stack with an options file that contains at least the following entries:

```
CMAKE_FLAGS+="-DDUNE_PYTHON_VIRTUALENV_SETUP=1
              -DDUNE_PYTHON_ALLOW_GET_PIP=1"
```

# Mathematical Abstractions

dune-structures provides solvers for elastostatic and elastodynamic simulations with a variety of linear and nonlinear material laws.

## Basic equations

Given a reference configuration $`T\subseteq\mathbb{R}^3`$ of an elastic body, we aim to calculate the displacement field $`\mathbf{u}:T\rightarrow\mathbb{R}^3`$ that describes the deformation of the body subject to the given forces.

## Strain Measures

Material laws in hyperelasticity are often expressed using a variety of different strain measures.
These are summarized here with the goal to document the implemented physics:

* *Deformation gradient*: $`\mathbf{F} = \mathbb{I} + \nabla\mathbf{u}`$
* *Infinitesimal strain*: $`\frac{1}{2}(\nabla\mathbf{u} + \nabla\mathbf{u}^T)`$
* *Right Cauchy-Green tensor*: $`\mathbf{C} = \mathbf{F}^T\mathbf{F}`$
* *Left Cauchy-Green tensor*: $`\mathbf{B} = \mathbf{F}\mathbf{F}^T`$
* *Cauchy-Green strain tensor*: $`\frac{1}{2}(\mathbf{C}-\mathbb{I})`$
* The *Cauchy-Green invariants*:
	* $`I_1 = tr(\mathbf{B})`$
	* $`I_2 = \frac{1}{2}(tr(\mathbf{B})^2 - tr(\mathbf{B}\mathbf{B})`$
	* $`I_3 = \det \mathbf{B}`$
* The *isochoric Cauchy-Green invariants*:
	* $`\bar{I}_1 = J^{-\frac{2}{3}}I_1`$
	* $`\bar{I}_2 = J^{-\frac{4}{3}}I_2`$
	* $`J=\sqrt{I_3}`$

## Finite Element formulation

# Software Design

# Solver Component Documentation

# How to...

here are some guidelines for likely extensions of the dune-structures code base.

## ... add a material law

## ... add a solver step

## ... do XYZ

Write an email to [dominic.kempf@iwr.uni-heidelberg.de](mailto:dominic.kempf@iwr.uni-heidelberg.de).