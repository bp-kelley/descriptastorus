from __future__ import print_function
import unittest
from rdkit.Chem import AllChem
import numpy, math
from descriptastorus import make_store, DescriptaStore
from descriptastorus.descriptors import rdDescriptors, DescriptorGenerator
import contextlib, tempfile, os, shutil, sys
import datahook

one_smiles = "c1ccccc1 0"
many_smiles = "\n".join( [ "C"*i + "c1ccccc1 " + str(i) for i in range(10) ] )
from expected_chiral_data import expected_chiral_data

class TestCase(unittest.TestCase):
    def testChiralMorgans(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
            print("\n\nfilename:", fname, file=sys.stderr),
            print("storefilename:", storefname, file=sys.stderr),
            with open(fname, 'w') as f:
                f.write(many_smiles),
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="MorganChiral3Counts",
                                                  index_inchikey=True )
            make_store.make_store(opts)
            generator = DescriptorGenerator.REGISTRY["MorganChiral3Counts".lower()]
            with contextlib.closing(DescriptaStore(storefname)) as store:
                for i in range(10):
                    r = store.descriptors().get(i)
                    self.assertEqual(r, expected_chiral_data[i])
                    
                
        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)
                
    def testChiralMorgans2(self):
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
                                                  seperator=" ", descriptors="MorganChiral3Counts",
                                                  index_inchikey=False )
            make_store.make_store(opts)

            with contextlib.closing(DescriptaStore(storefname)) as store:
                for i in range(10):
                    r = store.descriptors().get(i)
                    self.assertEqual(r, expected_chiral_data[i])
                
        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)


