This module contains the source code for the simulation of single cell mechanics and stress fibre formation executed in the STRUCTURES cluster of excellence at Heidelberg university. For questions write to [dominic.kempf@iwr.uni-heidelberg.de](mailto:dominic.kempf@iwr.uni-heidelberg.de).

# Installing the software

There is two principal ways to use this software: As dockerized applications or by building it from source. The former is recommended if you intend to *run* simulations, the latter is necessary if you want to *develop* your own simulators.

## Dockerized setup

The dockerized setup runs with only minimal requirements on the host machine:

* Docker (e.g. installed by following [these instructions from Docker](https://docs.docker.com/install/linux/docker-ce/debian/#install-using-the-repository)).
* docker-compose (available from the same repository)

These two commands will spin up the dune-structures frontend:

```
wget https://gitlab.dune-project.org/dominic/dune-structures/-/raw/master/docker-compose.yml
docker-compose up
```

You can now head to `0.0.0.0:5000` in your browser and use the visual programming frontend.

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
* [dune-blocklab](https://gitlab.dune-project.org/dominic/dune-blocklab.git)

Note that `dune-codegen` requires special care in setup. It needs to
* be cloned with the `--recursive` flag
* You need to run `./patches/apply_patches.sh` from the module root directory after cloning

Similarly, `dune-blocklab` and `dune-structures` need to be cloned recursively.

Currently, the latest versions (master) of all Dune modules are necessary.
Furthermore, dune-structures requires the following external packages:

* A Python interpreter >= 3.6 (e.g. Debian packages `python3-dev`)
* pip (e.g. the Debian package `python3-pip`)
* MuParser (e.g. the Debian package `libmuparser-dev`)
* Gmsh (e.g. Debian packages `gmsh`)
* `yaml-cpp` e.g. by installing the Debian package `libyaml-cpp-dev`
* The Python packages `cerberus` and `pyaml`

You should build your Dune stack with an options file that contains at least the following entries:

```
CMAKE_FLAGS+="-DDUNE_PYTHON_VIRTUALENV_SETUP=1
              -DDUNE_PYTHON_ALLOW_GET_PIP=1"
```

# How to...

here are some guidelines for likely extensions of the dune-structures code base.

## ... add a material law

- Have a look at the file `python/dune/structures/material.py`.
- Add a class to that file that inherits from `MaterialLawBase` and implements its fields
- Add the ID and material names you added in python to the `std::map` data structure in `dune/structures/material.hh`
- Register your material for code generation by adding it to the list of materials in `operators/operators.ufl`

## ... add a solver block

This documentation has moved entirely to dune-blocklab.

## ... add a new application

Just have a look at the applications in the `apps` subdirectory and add one that has the
grid implementation, finite element implementation you want and that includes all the blocks
you want. Watch out for `apps/CMakeLists.txt` to see how you integrate your app into the frontend.

## ... do XYZ

Write an email to [dominic.kempf@iwr.uni-heidelberg.de](mailto:dominic.kempf@iwr.uni-heidelberg.de).
