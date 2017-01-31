#!/usr/bin/env python
"""storus
see Description below
"""
from __future__ import print_function
from descriptastorus import DescriptaStore
from descriptastorus.descriptors import MakeGenerator
import argparse, os, shutil, time, random

import sys
from rdkit import rdBase
rdBase.DisableLog("rdApp.*")


parser = argparse.ArgumentParser()
parser.add_argument("storage",
                    help="directory in which to store the descriptors")

parser.add_argument("--descriptors", default="Morgan3Counts,RDKit2D")
parser.add_argument("--samples", default=1000, type=int)

opts = parser.parse_args()
store = DescriptaStore(opts.storage)

N = len(store)
gen = MakeGenerator(opts.descriptors.split(","))
import math

for i in range(opts.samples):
    idx = random.randint(0,N)
    
    v = store.descriptors().get(idx)
    smiles = store.molIndex().getMol(idx)

    res = gen.process(smiles)
    if res is None:
        assert v == tuple([0]*len(v))
        assert smiles == "NOSTRUCT"
        continue
    
    data = []
    
    for x in gen.process(smiles):
        if math.isnan(x): data.append(None)
        else: data.append(x)
    v2 = []
    for x in v:
        if math.isnan(x): v2.append(None)
        else: v2.append(x)
    assert v2 == data, "idx:%s smiles:%s\n%r\n\t%r"%(idx, smiles, v, data)

        
            
            

