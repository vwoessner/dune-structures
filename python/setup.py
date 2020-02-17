#!/usr/bin/env python

from setuptools import setup

setup(name='dune.structures',
      version='0.1',
      namespace_packages=['dune'],
      description='Cell mechanics simulation tools for the STRUCTURES project',
      author='Dominic Kempf <dominic.kempf@iwr.uni-heidelberg.de>',
      url='https://gitlab.dune-project.org/dominic/dune-structures.git',
      packages=['dune.structures',
                ],
      install_requires=['dune.codegen',
                        'pyaml',
                        'pygmsh',
                        'meshio',
                        ],
      entry_points = {
        "console_scripts": [
            "generate_cell_mesh = dune.structures.gmsh:entrypoint_generate_mesh",
            "generate_tangential_derivatives = dune.structures.diffgeo:generate_tangential_derivatives",
        ]
    })

