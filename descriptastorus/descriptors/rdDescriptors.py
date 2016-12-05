"""
This is a canned example for generating descriptors
Please modify as necessary
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rd
import numpy
from .DescriptorGenerator import DescriptorGenerator

import sys

def clip(v, name):
    if v > 255:
        print >> sys.stderr, "%s has >255 counts for a feature vector"%name
        v = 255
    return v
        

class RDKitMorgan3CountsAndDescriptors(DescriptorGenerator):
    def __init__(self, props=rd.Properties()):
        # specify names and numpy types for all columns
        morgan = [("m3-%d"%d, numpy.uint8) for d in range(2048)]
        descriptors = [(p, numpy.float64) for p in props.GetPropertyNames()]
        self.columns = morgan + descriptors
        self.props = props

    def GetColumns(self):
        """Returns [(name, numpy.dtype), ...] for all columns being computed"""
        return self.columns

    def processMol(self, m, smiles):
        counts = list(rd.GetHashedMorganFingerprint(m, radius=3, nBits=2048))
        counts = [ clip(x,smiles) for x in counts ]
        counts.extend(self.props.ComputeProperties(m))
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

