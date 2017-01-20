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
        print >> sys.stderr, "%s has >255 counts for a feature vector"%name
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
        self.columns = morgan

    def GetColumns(self):
        """Returns [(name, numpy.dtype), ...] for all columns being computed"""
        return self.columns

    def processMol(self, m, smiles):
        counts = list(rd.GetHashedMorganFingerprint(m,
                                                    radius=self.radius, nBits=self.nbits))
        counts = [ clip(x,smiles) for x in counts ]
        return counts
        
    def process(self, smiles):
        """smiles
        generate descriptors from a smiles string using the specified
        properties.  

        Default is to return morgan3 folded counts clipped to 255 and
        use rdkit 2D properties.
        """
        try:
            m = Chem.MolFromSmiles(smiles)
        except:
            return None

        if m == None:
            return None

        return self.processMol(m, smiles)

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
            columns = self.columns = []
            for p in properties:
                if p in sorted(FUNCS):
                    columns.append((name, numpy.float64))
                else:
                    raise KeyError("Unable to find specified property %s"%p)
        
    def GetColumns(self):
        """Returns [(name, numpy.dtype), ...] for all columns being computed"""
        return self.columns

    def processMol(self, m, smiles):
        return [ applyFunc(name, m) for name, _ in self.columns ]
    

RDKit2D()
