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
        logging.warning("Extra columns in store1:\n\t%s", "\n\t".join(s1-s2))
        logging.warning("Extra columns in store2:\n\t%s", "\n\t".join(s2-s1))
        
    for i,(a,b) in enumerate(zip(store1.getDescriptorNames(),
                                 store2.getDescriptorNames())):
        if a != b:
            logging.warning("First differing element at index %s (%s != %s)", i,a,b)
            break
        
    sys.exit(1)
    
# check stores...    

