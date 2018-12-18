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
#!/usr/bin/env python
"""storus
see Description below
"""
from __future__ import print_function
from descriptastorus import DescriptaStore
from descriptastorus.descriptors import MakeGenerator
import argparse, os, shutil, time, random

import math, sys
from rdkit import rdBase
rdBase.DisableLog("rdApp.*")


parser = argparse.ArgumentParser()
parser.add_argument("storage",
                    help="directory in which to store the descriptors")

parser.add_argument("--samples", default=1000, type=int)

def main():
    opts = parser.parse_args()
    store = DescriptaStore(opts.storage)
    
    N = len(store)
    gen = store.getDescriptorCalculator()
    
    
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
            idx = random.randint(0,N-1)
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
    
            
                
                
    
    
    
