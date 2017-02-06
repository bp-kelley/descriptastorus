"""
This is a canned example for generating descriptors
Please modify as necessary
"""
from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors as rd
import numpy
from .DescriptorGenerator import DescriptorGenerator
import logging

import sys

def clip(v, name):
    if v > 255:
        logging.warning("%s has >255 counts for a feature vector clipping", name)
        v = 255
    return v
        

class MorganCounts(DescriptorGenerator):
    """Computes Morgan3 bitvector counts"""
    NAME = "Morgan%sCounts"
    def __init__(self, radius=3, nbits=2048):
        if radius == 3 and nbits == 2048:
            self.NAME = self.NAME % "3"
        else:
            self.NAME = self.NAME % ("%s-%s"%radius,nbits)
            
        DescriptorGenerator.__init__(self)
        # specify names and numpy types for all columns
        self.radius = radius
        self.nbits = nbits
        morgan = [("m3-%d"%d, numpy.uint8) for d in range(nbits)]
        self.columns += morgan

    def calculateMol(self, m, smiles, internalParsing=False):
        counts = list(rd.GetHashedMorganFingerprint(m,
                                                    radius=self.radius, nBits=self.nbits))
        counts = [ clip(x,smiles) for x in counts ]
        return counts        

MorganCounts()

FUNCS = {name:func for name, func in Descriptors.descList}
def applyFunc(name, m):
    try:
        return FUNCS[name](m)
    except:
        logging.exception("function application failed (%s->%s)",
            name, Chem.MolToSmiles(m))
                       
        return 0.0

            
class RDKit2D(DescriptorGenerator):
    """Computes all RDKit Descriptors"""
    NAME = "RDKit2D"
    def __init__(self, properties=None):
        DescriptorGenerator.__init__(self)
        # specify names and numpy types for all columns
        if not properties:
            self.columns = [ (name, numpy.float64) for name,func in sorted(Descriptors.descList) ]
        else:
            columns = self.columns
            for name in properties:
                if name in sorted(FUNCS):
                    columns.append((name, numpy.float64))
                else:
                    raise KeyError("Unable to find specified property %s"%name)
        
    def calculateMol(self, m, smiles, internalParsing=False):
        res = [ applyFunc(name, m) for name, _ in self.columns ]
        return res
    

RDKit2D()

    
