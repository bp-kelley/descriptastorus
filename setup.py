#!/usr/bin/env python
#  Copyright (c) 2018, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
import sys
try:
  import rdkit
except:
  print("Descriptastorus requires rkdit to function, this is not installable by pip", file=sys.stderr)
  print(" see https://rdkit.org for more information", file=sys.stderr)
  sys.exit(1)
  
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
  VERSION="2.2.0" # hardcode version

setup(name='descriptastorus',
      version=VERSION,
      description='Descriptor storage and molecular file indexing',
      author='Brian Kelley',
      author_email='brian.kelley@novartis.com',
      url='https://bitbucket.org/novartisnibr/rdkit-descriptastorus/',
      install_requires=['pandas_flavor'],
      test_suite='nose.collector',
      tests_require=['nose', 'pandas_flavor'],
      include_package_data=True,
      entry_points={
          'console_scripts': [
              'storus = descriptastorus.cli.storus:main',
              'storus-validate = descriptastorus.cli.validate:main',
              ]
          },
                                                            
      packages = find_packages())

