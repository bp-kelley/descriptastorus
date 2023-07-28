from __future__ import print_function
import unittest
from rdkit.Chem import AllChem
import numpy, math
from descriptastorus import make_store, DescriptaStore
from descriptastorus.descriptors.DescriptorGenerator import DescriptorGenerator
import contextlib, tempfile, os, shutil, sys
import datahook
make_store.DEFAULT_KEYSTORE = "dbmstore"

one_smiles = "c1ccccc1 0"
many_smiles = "\n".join( [ "C"*i + "c1ccccc1 " + str(i) for i in range(10) ] )

# this shoudl throw type errors in MakeStore
class NanDescriptors(DescriptorGenerator):
    NAME="NANDescriptors"
    def GetColumns(self):
        return [('a', numpy.int32),
                ('b', numpy.float32),
                ('c', numpy.float64),
                ('d', numpy.uint8)]
    
    def processMol(self, m, smiles, internalParsing=False):
        return [None]*4
    
NanDescriptors()

# this should store a calculated flag error
class NanDescriptorsWithCalcFlags(DescriptorGenerator):
    NAME="NANDescriptorsWithCalcFlags"
    def __init__(self):
        DescriptorGenerator.__init__(self)
        self.columns =[('a', numpy.int32),
                       ('b', numpy.float32),
                       ('c', numpy.float64),
                       ('d', numpy.uint8)]
    
    def calculateMol(self, m, smiles, internalParsing):
        return [None]*4
    
descriptors = NanDescriptorsWithCalcFlags()

class TestCase(unittest.TestCase):
    def testRawNones(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
            with open(fname, 'w') as f:
                f.write(one_smiles)
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="NANDescriptors",
                                                  index_inchikey=True )
            try:
                make_store.make_store(opts)
                self.assertFalse(True) # should be a type error
            except TypeError as e:
                self.assertTrue("For column" in str(e))
                self.assertTrue("can't convert" in str(e))
                
        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)

    def testDescriptorWithNones(self):
        res = descriptors.processMol(None, "", False)
        self.assertEqual(res, [False, 0, 0, 0.0, 0.0])
        self.assertFalse(res[0])
        self.assertFalse(None in res)
        
        
    def testNonesWithCalcFlags(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
            with open(fname, 'w') as f:
                f.write(one_smiles)
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="NANDescriptorsWithCalcFlags",
                                                  index_inchikey=True )
            make_store.make_store(opts)
            with contextlib.closing(DescriptaStore(storefname)) as store:
                self.assertFalse( store.descriptors().get(0)[0] )
                
        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)

if __name__ == '__main__':  #pragma: no cover
    unittest.main()
