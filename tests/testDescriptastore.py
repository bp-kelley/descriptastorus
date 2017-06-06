from __future__ import print_function
import unittest
from rdkit.Chem import AllChem
from descriptastorus import make_store, DescriptaStore
from descriptastorus.descriptors.rdDescriptors import RDKit2D

import contextlib, tempfile, os, shutil, sys
import datahook

one_smiles = "c1ccccc1 0"
many_smiles = "\n".join( [ "C"*i + "c1ccccc1 " + str(i) for i in range(10) ] + ["NOSTRUCT foo"] )

from rdkit.Chem import Descriptors

class RDKit2DSubset(RDKit2D):
    NAME="RDKit2DSubset"
    def __init__(self):
        RDKit2D.__init__(self, properties=[
            'ExactMolWt',
            'NumAliphaticRings', 'NumAromaticCarbocycles',
            'NumAromaticHeterocycles', 'NumAromaticRings'])
RDKit2DSubset()

def toDict( v ):
    return {n:x for n,x in zip([
            'RDKit2DSubset_calculated', 'ExactMolWt',
            'NumAliphaticRings', 'NumAromaticCarbocycles',
            'NumAromaticHeterocycles', 'NumAromaticRings'],
                               v)}

class TestCase(unittest.TestCase):
    def testOffByOne(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
            print("\n\nfilename:", fname, file=sys.stderr)
            print("storefilename:", storefname, file=sys.stderr)
            with open(fname, 'w') as f:
                f.write(one_smiles)
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  batchsize=1,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="RDKit2DSubset",
                                                  index_inchikey=True )
            make_store.make_store(opts)

            with contextlib.closing(DescriptaStore(storefname)) as store:
                self.assertEqual( store.lookupName("0"), 0)

                self.assertEqual( store.lookupInchiKey("UHOVQNZJYSORNB-UHFFFAOYSA-N"), [0])
                self.assertEqual(store.descriptors().get(0), (True, 78.046950192, 0.0, 1.0, 0.0, 1.0))

            try:
                store.lookupInchiKey("MY DOG HAS FLEAS")
                self.assertTrue(False) # should not get here
            except KeyError:
                pass

        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)
                
    def testMany(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
            with open(fname, 'w') as f:
                f.write(many_smiles)
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="RDKit2DSubset",
                                                  index_inchikey=True )
            make_store.make_store(opts)

            with contextlib.closing(DescriptaStore(storefname)) as store:

                for i in range(10):
                    self.assertEqual( store.lookupName(str(i)), i)

                self.assertEqual(store.descriptors().get(0), (True, 78.046950192, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(1), (True, 92.062600256, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(2), (True, 106.07825032, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(3), (True, 120.093900384, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(4), (True, 134.109550448, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(5), (True, 148.125200512, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(6), (True, 162.140850576, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(7), (True, 176.15650064, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(8), (True, 190.172150704, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(9), (True, 204.187800768, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(10), (False, 0.0, 0.0, 0.0, 0.0, 0.0))

                self.assertEqual(store.descriptors().getDict(7), toDict((True, 176.15650064, 0.0, 1.0, 0.0, 1.0)))

                calc = store.getDescriptorCalculator()
                
                for i in range(10):
                    m = store.molIndex().getRDMol(i)
                    sm = AllChem.MolToSmiles(m)
                    inchi = AllChem.InchiToInchiKey(AllChem.MolToInchi(m))
                    self.assertEqual( store.lookupInchiKey(inchi), [i])
                    v = store.descriptors().get(i)
                    sv = tuple(calc.process(sm))
                    self.assertEqual(v, sv)

        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)

    def testManyNoInchi(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
            with open(fname, 'w') as f:
                f.write(many_smiles)
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  batchsize=1,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="RDKit2DSubset",
                                                  index_inchikey=False )
            make_store.make_store(opts)

            origdata = many_smiles.split("\n")
            
            with contextlib.closing(DescriptaStore(storefname)) as store:
                for i in range(10):
                    self.assertEqual( store.lookupName(str(i)), i)

                self.assertEqual(store.descriptors().get(0), (True, 78.046950192, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(1), (True, 92.062600256, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(2), (True, 106.07825032, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(3), (True, 120.093900384, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(4), (True, 134.109550448, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(5), (True, 148.125200512, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(6), (True, 162.140850576, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(7), (True, 176.15650064, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(8), (True, 190.172150704, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(9), (True, 204.187800768, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual(store.descriptors().get(10), (False, 0.0, 0.0, 0.0, 0.0, 0.0))                
                for i in range(10):
                    m = store.molIndex().getRDMol(i)
                    smiles, name = store.molIndex().get(i)
                    self.assertEqual(name, str(i))
                    self.assertEqual(smiles, origdata[i].split()[0])


        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)
                

    def testContainer(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
            with open(fname, 'w') as f:
                f.write(many_smiles)
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="RDKit2DSubset,RDKit2DSubset",
                                                  index_inchikey=True )
            make_store.make_store(opts)

            with contextlib.closing(DescriptaStore(storefname)) as store:

                for i in range(10):
                    self.assertEqual( store.lookupName(str(i)), i)

                for i in range(10):
                    m = store.molIndex().getRDMol(i)
                    inchi = AllChem.InchiToInchiKey(AllChem.MolToInchi(m))
                    self.assertEqual( store.lookupInchiKey(inchi), [i])
                self.assertEqual(store.descriptors().get(0), (True, 78.046950192, 0.0, 1.0, 0.0, 1.0)*2)
                self.assertEqual(store.descriptors().get(1), (True, 92.062600256, 0.0, 1.0, 0.0, 1.0)*2)
                self.assertEqual(store.descriptors().get(2), (True, 106.07825032, 0.0, 1.0, 0.0, 1.0)*2)
                self.assertEqual(store.descriptors().get(3), (True, 120.093900384, 0.0, 1.0, 0.0, 1.0)*2)
                self.assertEqual(store.descriptors().get(4), (True, 134.109550448, 0.0, 1.0, 0.0, 1.0)*2)
                self.assertEqual(store.descriptors().get(5), (True, 148.125200512, 0.0, 1.0, 0.0, 1.0)*2)
                self.assertEqual(store.descriptors().get(6), (True, 162.140850576, 0.0, 1.0, 0.0, 1.0)*2)
                self.assertEqual(store.descriptors().get(7), (True, 176.15650064, 0.0, 1.0, 0.0, 1.0)*2)
                self.assertEqual(store.descriptors().get(8), (True, 190.172150704, 0.0, 1.0, 0.0, 1.0)*2)
                self.assertEqual(store.descriptors().get(9), (True, 204.187800768, 0.0, 1.0, 0.0, 1.0)*2)
                self.assertEqual(store.descriptors().get(10), (False, 0.0, 0.0, 0.0, 0.0, 0.0)*2)                

        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)
                
