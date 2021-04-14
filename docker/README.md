# Docker Images

This directory contains source files for building Docker images required to
build and/or test the `dune-structures` module.

## Building the environment image

The environment image contains all packages for building and testing `dune-structures`.
It must be built manually using the `build-env.dockerfile`.
To use it in the GitLab CI/CD, it must be pushed to the container registry.
Follow these steps:

```bash
cd docker
docker build -f build-env.dockerfile -t registry.dune-project.org/lukas.riedel/dune-structures:env .
docker login registry.dune-project.org
docker push registry.dune-project.org/lukas.riedel/dune-structures:env
```
