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

parser.add_argument("--samples", default=1000, type=int)

opts = parser.parse_args()
store = DescriptaStore(opts.storage)

N = len(store)
gen = MakeGenerator(store.options['descriptors'].split(","))
import math

randomize=True
if opts.samples == -1:
    randomize=False
    opts.samples = len(store)

next = .05
for i in range(opts.samples):
    if i and float(i)/opts.samples > next:
        print("Validated %2.2f%%"%(next*100))
        next += .05
    if randomize:
        idx = random.randint(0,N)
    else:
        idx = i
    
    v = store.descriptors().get(idx)
    smiles = store.molIndex().getMol(idx)
    name = None
    try:
        name = store.molIndex().getName(idx)
    except:
        pass
    
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
    assert v2 == data, "idx:%s smiles:%s name:%s \n%r\n\t%r"%(idx, smiles, name,
                                                              v, data)

        
            
            

