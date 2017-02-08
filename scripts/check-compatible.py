#!/usr/bin/env python
"""storus
see Description below
"""
from __future__ import print_function
from descriptastorus import DescriptaStore
import argparse, os, logging, shutil, time, random

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
    logging.warning("Storages are not compatible, columns are different")
    s1 = set(store1.getDescriptorNames())
    s2 = set(store2.getDescriptorNames())
    
    if s1 == s2:
        logging.warning("Columns are the same but reordered")
    else:
        logging.warning("Extra columns in store1:\n%s", "\n\t".join(s1-s2))
        logging.warning("Extra columns in store2:\n%s", "\n\t".join(s2-s1))
        
    for i,(a,b) in enumerate(zip(store1.getDescriptorNames(),
                                 store2.getDescriptorNames())):
        if a != b:
            logging.warning("First differing element at index %s (%s != %s)", i,a,b)
            break
        
    sys.exit(1)
    
# check stores...    

