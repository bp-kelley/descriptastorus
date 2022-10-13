from __future__ import print_function
import unittest
from rdkit.Chem import AllChem
import numpy, math
from rdkit.DataStructs import IntSparseIntVect
from descriptastorus import make_store, DescriptaStore
from descriptastorus.descriptors import rdDescriptors, DescriptorGenerator
import contextlib, tempfile, os, shutil, sys
import datahook

one_smiles = "c1ccccc1 0"
many_smiles = "\n".join( [ "C"*i + "c1ccccc1 " + str(i) for i in range(10) ] )
from expected_FPdescriptor_data  import expected_chiral_data, expected_RDKFP_data, expected_AtomPair_data, expected_FeatureMorgan_data

#Some test smiles from the Huuskonen data set
testSmiles = {
'Natamycin': 'O1C(C)C(O)C(N)C(O)C1OC2C=CC=CC=CC=CCC(C)OC(=O)C=CC(O3)C3CC(O)CC4(O)OC(C(C(=O)O)C(O)C4)C2',
'n-propylbenzene': 'c1ccccc1CCC',
'D-limonene': 'C(=CCC(C(=C)C)C1)(C1)C',
'1-heptyne': 'CCCCCC#C'
}


class TestCase(unittest.TestCase):
    def testChiralMorgans(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
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
                    if r != expected_chiral_data[i]:
                        asdf
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
                    if r != expected_chiral_data[i]:
                        print(i)
                        print(r)
                        print("*"*33)
                        print(expected_chiral_data[i])
                        asdf

                    self.assertEqual(r, expected_chiral_data[i])
                
        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)
                
    def testRDKitFPBits(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
            with open(fname, 'w') as f:
                f.write("\n".join(['{0} {1}'.format(v,k) for k,v in testSmiles.items()]))
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="RDKitFPBits",
                                                  index_inchikey=False )
            make_store.make_store(opts)

            with contextlib.closing(DescriptaStore(storefname)) as store:
                for i in range(4):
                    r = store.descriptors().get(i)
                    self.assertEqual(r, expected_RDKFP_data[i])
                
        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)


    def testAtomPairCounts(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
            with open(fname, 'w') as f:
                f.write("\n".join(['{0} {1}'.format(v,k) for k,v in testSmiles.items()]))
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="AtomPairCounts",
                                                  index_inchikey=False )
            make_store.make_store(opts)

            with contextlib.closing(DescriptaStore(storefname)) as store:
                for i in range(4):
                    r = store.descriptors().get(i)
                    #print(r)
                    self.assertEqual(r, expected_AtomPair_data[i])
                
        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)

    def testFeatureMorganCounts(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
            with open(fname, 'w') as f:
                f.write("\n".join(['{0} {1}'.format(v,k) for k,v in testSmiles.items()]))
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="MorganFeature3Counts",
                                                  index_inchikey=False )
            make_store.make_store(opts)

            with contextlib.closing(DescriptaStore(storefname)) as store:
                for i in range(4):
                    r = store.descriptors().get(i)
                    if r !=  expected_FeatureMorgan_data[i]:
                        print(i)
                        print(r)
                        print(expected_FeatureMorgan_data[i])
                        asdf
                    #print(r)expected_FeatureMorgan_data[i]
                    self.assertEqual(r, expected_FeatureMorgan_data[i])
                
        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)

    def test_clip(self):
        nbits = 2048
        v = IntSparseIntVect(2048)
        for i in range(nbits):
            v[i] = i
            
        l = rdDescriptors.clip_sparse(v, 2048)
        for i,v in enumerate(l):
            assert v == min(i,255)

if __name__ == '__main__':  #pragma: no cover
    unittest.main()

