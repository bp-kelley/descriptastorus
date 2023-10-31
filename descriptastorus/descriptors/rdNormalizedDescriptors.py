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
from . import rdDescriptors
from . import dists
from . import hists
from collections import namedtuple
import scipy.stats as st
import numpy as np
import logging
from bisect import bisect

logger = logging.getLogger("descriptastorus")

cdfs = {}

for name, (dist, params, minV,maxV,avg,std) in dists.dists.items():
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    if dist in ['gilbrat', 'gibrat']:
        # fix change in scikit learn
        if hasattr(dist, 'gilbrat'):
            dist = 'gilbrat'
        else:
            dist = 'gibrat'
            
    dist = getattr(st, dist)
        
    
    # make the cdf with the parameters
    def cdf(v, dist=dist, arg=arg,loc=loc,scale=scale,minV=minV,maxV=maxV):
        v = dist.cdf(np.clip(v, minV, maxV), loc=loc, scale=scale, *arg)
        return np.clip(v, 0., 1.)
    
    cdfs[name] = cdf

#for name in rdDescriptors.FUNCS:
#    if name not in cdfs:
#        logger.warning("No normalization for %s", name)

def applyNormalizedFunc(name, m):
    if name not in cdfs:
        return 0.0
    try:
        return cdfs[name](rdDescriptors.applyFunc(name,m))
    except:
        logger.exception("Could not compute %s for molecule", name)
        return 0.0

class RDKit2DNormalized(rdDescriptors.RDKit2D):
    """Distribution normalized descriptors.
    These are then converted to a CDF so for a given value v
    the result is the percentage of the population under that value
    """
    NAME = "RDKit2DNormalized"

    def calculateMol(self, m, smiles, internalParsing=False):
        res = [ applyNormalizedFunc(name, m) for name, _ in self.columns ]
        return res   
    
RDKit2DNormalized()

histcdfs = {}
for name, bins in hists.hists.items():
    # make the cdf with the parameters
    def histcdf(v, bins=bins):
        p = bisect(bins, (v,))
        if p < len(bins):
            return bins[p][1]
        return 1.0
    
    histcdfs[name] = histcdf

def applyHistogramNormalizedFunc(name, m):
    if name not in cdfs:
        return 0.0
    try:
        return histcdfs[name](rdDescriptors.applyFunc(name,m))
    except:
        logging.exception("Could not compute %s for molecule", name)
        return 0.0

class RDKit2DHistogramNormalized(rdDescriptors.RDKit2D):
    """Histogrram normalized descriptors.
    These are then converted to a CDF so for a given value v
    the result is the percentage of the population under that value

    These are WAY faster than the corresponding Distribution normalized 
    descriptors
    """
    NAME = "RDKit2DHistogramNormalized"

    def calculateMol(self, m, smiles, internalParsing=False):
        res = [ applyHistogramNormalizedFunc(name, m) for name, _ in self.columns ]
        return res   
    
RDKit2DHistogramNormalized()
