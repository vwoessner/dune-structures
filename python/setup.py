#!/usr/bin/env python

from setuptools import setup

setup(name='dune.structures',
      version='0.1',
      namespace_packages=['dune'],
      description='Performance optimizing form compiler for the Dune project',
      author='Dominic Kempf <dominic.kempf@iwr.uni-heidelberg.de>',
      url='https://gitlab.dune-project.org/dominic/dune-structures.git',
      packages=['dune.structures',
                ],
      install_requires=['dune.codegen'],
      entry_points = {
        "console_scripts": [
        ]
    })

