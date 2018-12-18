from . import rdDescriptors
from . import dists
from collections import namedtuple
import scipy.stats as st
import numpy as np
import logging

cdfs = {}

for name, (dist, params, minV,maxV,avg,std) in dists.dists.items():
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    dist = getattr(st, dist)
    
    # make the cdf with the parameters
    def cdf(v, dist=dist, arg=arg,loc=loc,scale=scale,minV=minV,maxV=maxV):
        v = dist.cdf(np.clip(v, minV, maxV), loc=loc, scale=scale, *arg)
        return np.clip(v, 0., 1.)
    
    cdfs[name] = cdf

for name in rdDescriptors.FUNCS:
    if name not in cdfs:
        logging.warning("No normalization for %s", name)

def applyNormalizedFunc(name, m):
    if name not in cdfs:
        return 0.0
    return cdfs[name](rdDescriptors.applyFunc(name,m))

class RDKit2DNormalized(rdDescriptors.RDKit2D):
    NAME = "RDKit2DNormalized"

    def calculateMol(self, m, smiles, internalParsing=False):
        res = [ applyNormalizedFunc(name, m) for name, _ in self.columns ]
        return res   
    
RDKit2DNormalized()
