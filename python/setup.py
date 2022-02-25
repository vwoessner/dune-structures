#!/usr/bin/env python

from setuptools import setup

setup(
    name="dune.structures",
    version="0.1",
    namespace_packages=["dune"],
    description="Cell mechanics simulation tools for the STRUCTURES project",
    author="Dominic Kempf <dominic.kempf@iwr.uni-heidelberg.de>",
    url="https://gitlab.dune-project.org/dominic/dune-structures.git",
    packages=[
        "dune.structures",
    ],
    install_requires=[
        "dune.codegen",
        "pyaml",
        "pygmsh",
        "meshio",
        "vtk",
        "numpy",
        "pytest",
        "ruamel.yaml",
        "scipy",
        "pandas",
        "seaborn",
        "meshio",
    ],
    entry_points={
        "console_scripts": [
            "generate_cell_mesh = dune.structures.gmsh:entrypoint_generate_mesh",
            "generate_tangential_derivatives = dune.structures.diffgeo:generate_tangential_derivatives",
            "structures = dune.structures.cli:cli",
            "plot_histogram = dune.structures.plotting.histograms:entrypoint",
            "plot_population = dune.structures.plotting.population:entrypoint",
            "plot_mean = dune.structures.plotting.displacement:entrypoint_mean",
            "plot_displacement = dune.structures.plotting.displacement:entrypoint_single",
            "plot_density = dune.structures.plotting.density:entrypoint",
        ]
    },
)
