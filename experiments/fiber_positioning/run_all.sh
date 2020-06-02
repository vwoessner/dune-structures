#!/bin/bash

# Run the entire experiment

set -e

# Penalty factor in one interval
./run_experiments.py hansbo_beam.ini 10 0.5 0.6 100 1.0
./run_experiments.py hansbo_beam.ini 10 0.5 0.6 100 10.0
./run_experiments.py hansbo_beam.ini 10 0.5 0.6 100 100.0
./run_experiments.py hansbo_beam.ini 10 0.5 0.6 100 1000.0
./run_experiments.py hansbo_beam.ini 10 0.5 0.6 100 10000.0

# The full curve to see that it works fine all over the domain
./run_experiments.py hansbo_beam.ini 10 0.0 1.0 500 1000.0

# See impact of refinement on the profile
./run_experiments.py hansbo_beam.ini 1 0.0 1.0 100 1000.0
./run_experiments.py hansbo_beam.ini 5 0.0 1.0 100 1000.0
./run_experiments.py hansbo_beam.ini 20 0.0 1.0 100 1000.0