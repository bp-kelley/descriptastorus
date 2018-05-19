#!/bin/bash

set -e

python setup.py install --single-version-externally-managed --record=record.txt
nosetests -v
