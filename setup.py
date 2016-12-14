#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='descriptastorus',
      version='1.0',
      description='Descriptor storage and molecular file indexing',
      author='Brian Kelley',
      author_email='brian.kelley@novartis.com',
      url='https://bitbucket.org/novartisnibr/rdkit-descriptastorus/',
      test_suite='nose.collector',
      tests_require=['nose'],
      include_package_data=True,
      packages = find_packages())

