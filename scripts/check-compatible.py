#!/usr/bin/env python
"""storus
see Description below
"""
from __future__ import print_function
from descriptastorus import DescriptaStore
import argparse, os, shutil, time, random

import sys
from rdkit import rdBase
rdBase.DisableLog("rdApp.*")


parser = argparse.ArgumentParser()
parser.add_argument("storage",
                    help="directory in which to store the descriptors")
parser.add_argument("storage2",
                    help="directory in which to store the descriptors")

opts = parser.parse_args()
store1 = DescriptaStore(opts.storage)
store2 = DescriptaStore(opts.storage2)

if store1.getDescriptorNames() != store2.getDescriptorNames():
    logging.warning("Storages are not compatable, columns are different")
    sys.exit(1)
    
# check stores...    

