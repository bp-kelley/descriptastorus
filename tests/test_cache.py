from __future__ import print_function
import unittest
from rdkit.Chem import AllChem
import numpy, math
from descriptastorus import make_store, DescriptaStore
from descriptastorus.descriptors import (
    rdDescriptors, rdNormalizedDescriptors, DescriptorGenerator, dists)
import contextlib, tempfile, os, shutil, sys
import datahook
from descriptastorus.descriptors import DescriptorGenerator

class FakeD(DescriptorGenerator):
    NAME = "FAKED"
    def GetColumns(self):
        return ["FAKED"]
    
    def calculateMol(self, m, smiles, internalParsing):
        return [len(smiles)]

d = FakeD()


class TestCase(unittest.TestCase):
    def testCache(self):
        smiles = ["C" * i for i in range(1,40)]
        smiles2 = ["C" * i for i in range(40,90)]
        
        m1, res = d.processSmiles(smiles)
        self.assertEqual(len(d.cache), len(smiles))

        # all found
        m2, res2 = d.processSmiles(smiles)
        self.assertEqual(res, res2)
        self.assertEqual(d.cache_hit, len(smiles))
        self.assertEqual(m1, m2) # mols should be the same
        
        # some found
        d.cache_hit = d.cache_miss = 0
        _, res3 = d.processSmiles(smiles + smiles2)
        self.assertEqual(d.cache_hit, len(smiles))
        self.assertEqual(d.cache_miss, len(smiles2))

        d.cache.clear()
        d.cache_hit = d.cache_miss = 0
        _, res4 = d.processSmiles(smiles + smiles2)
        self.assertEqual(d.cache_hit, 0)
        self.assertEqual(d.cache_miss, len(smiles) + len(smiles2))

        self.assertEqual(res3, res4)

        smiles = ["X"] + ["C" * i for i in range(1,40)]
        d.cache_hit = d.cache_miss = 0

        m1, res = d.processSmiles(smiles)
        self.assertEqual(d.cache_miss, 1)
        self.assertEqual(res[0], None)
                
if __name__ == '__main__':  #pragma: no cover
    unittest.main()

