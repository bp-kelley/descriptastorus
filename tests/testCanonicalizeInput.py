from __future__ import print_function
import unittest
from rdkit.Chem import AllChem
import numpy, math
from descriptastorus import make_store, DescriptaStore
from descriptastorus.descriptors.DescriptorGenerator import DescriptorGenerator
import contextlib, tempfile, os, shutil, sys
import datahook

one_smiles = "c1ccccc1 0"
many_smiles = "\n".join( [ "C"*i + "c1ccccc1 " + str(i) for i in range(10) ] )

class Canonicalize(DescriptorGenerator):
    NAME="Canonicalize"
    canonicalCount = 0
    def GetColumns(self):
        return [('count', numpy.float32)]

    def molFromSmiles(self, smiles):
        self.canonicalCount += 1
        m = AllChem.MolFromSmiles(smiles)
        m.SetProp("ccount", str(self.canonicalCount))
        return m
    
    def processMol(self, m, smiles, internalParsing):
        assert internalParsing == True
        return [int(m.GetProp("ccount"))]
        
c = Canonicalize()

class TestCase(unittest.TestCase):
    def testCanonicalSmiles(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
            print("\n\nfilename:", fname, file=sys.stderr)
            print("storefilename:", storefname, file=sys.stderr)
            with open(fname, 'w') as f:
                f.write(many_smiles)
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="Canonicalize",
                                                  index_inchikey=True )
            make_store.make_store(opts)
            
            with contextlib.closing(DescriptaStore(storefname)) as store:
                counts = []
                for i in range(10):
                    r = store.descriptors().get(i)
                    self.assertEqual(r[0], i+1)


                
        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)
                
    def testCanonicalSmiles2(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
            print("\n\nfilename:", fname, file=sys.stderr)
            print("storefilename:", storefname, file=sys.stderr)
            with open(fname, 'w') as f:
                f.write(many_smiles)
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="Canonicalize",
                                                  index_inchikey=False )
            make_store.make_store(opts)
            
            with contextlib.closing(DescriptaStore(storefname)) as store:
                counts = []
                for i in range(10):
                    r = store.descriptors().get(i)
                    self.assertEqual(r[0], i+1)


                
        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)
