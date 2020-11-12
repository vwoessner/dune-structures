FROM registry.dune-project.org/dominic/dune-blocklab/frontend:latest

RUN duneci-install-module --recursive https://gitlab.dune-project.org/dominic/dune-structures.git

# Build the applications and their frontend specification
WORKDIR /duneci/modules/dune-structures/build-cmake
RUN make frontend

# Start gunicorn
CMD ["/duneci/modules/dune-structures/build-cmake/frontend.sh"]
