#!/usr/bin/env python

from setuptools import setup, find_packages
try:
  from commands import getstatusoutput
except ImportError:
  from subprocess import getstatusoutput

import logging
  
status, output = getstatusoutput("git describe --tags")

if not status:
    try:
        data = output.split("-")
        if len(data) == 1:
            VERSION = data[0]
        elif data and data[0].lower() == "release":
            VERSION = "%s.%s"%(data[1], data[2])
        else:
            VERSION = "%s.%s"%(data[0], data[1])
    except:
        raise RunTimeError("git tags must be in the form release-x.y.z or simply x.y.z")
else:
  raise RuntimeError("git must be accessible in your path and at least one release tag should exist in the repo, aborting")

setup(name='descriptastorus',
      version=VERSION,
      description='Descriptor storage and molecular file indexing',
      author='Brian Kelley',
      author_email='brian.kelley@novartis.com',
      url='https://bitbucket.org/novartisnibr/rdkit-descriptastorus/',
      test_suite='nose.collector',
      tests_require=['nose'],
      include_package_data=True,
      packages = find_packages())

