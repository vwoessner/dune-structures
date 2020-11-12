#!/bin/bash

#
# Build Docker images and push them to the registry
#

set -e

docker login registry.dune-project.org

# Build and push the 'frontend' image
image=registry.dune-project.org/dominic/dune-structures/frontend
docker build "$@" -t $image - < frontend.dockerfile
docker push $image

docker logout registry.dune-project.org
