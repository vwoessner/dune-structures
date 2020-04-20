#!/bin/bash

set -e

# The only argument to this script is the tag that should be applied to the image
tag="${1:latest}"
echo "Building docker container with tag $tag"

docker login registry.dune-project.org/dominic/dune-structures
docker build -t registry.dune-project.org/dominic/dune-structures:$tag
docker push registry.dune-project.org/dominic/dune-structures:$tag
docker logout registry.dune-project.org/dominic/dune-structures
