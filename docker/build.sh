#!/bin/bash

#
# Build the Docker image and push it to the registry
#

set -e


docker login registry.dune-project.org

for toolchain in gcc-9-17 clang-9-libcpp-17
do
  image=registry.dune-project.org/dominic/dune-structures/base:debian-11-$toolchain
  docker build --pull "$@" --build-arg TOOLCHAIN=$toolchain -t $image - < base.dockerfile
  docker push $image
done

# Build and push the 'frontend' image
image=registry.dune-project.org/dominic/dune-structures/frontend
docker build "$@" -t $image - < frontend.dockerfile
docker push $image

docker logout registry.dune-project.org
