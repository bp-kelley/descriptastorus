from __future__ import print_function
import unittest
from rdkit.Chem import AllChem
from descriptastorus import make_store, DescriptaStore
import tempfile, os, shutil, sys
import datahook

one_smiles = "c1ccccc1 0"
many_smiles = "\n".join( [ "C"*i + "c1ccccc1 " + str(i) for i in range(10) ] )

class TestCase(unittest.TestCase):
    def testOffByOne(self):
        try:
            fname = tempfile.mktemp()
            storefname = tempfile.mktemp()
            with open(fname, 'w') as f:
                f.write(one_smiles)
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="RDKit2D",
                                                  index_inchikey=True )
            make_store.make_store(opts)

            store = DescriptaStore(storefname)
                
            self.assertEqual( store.lookupName("0"), 0)

            self.assertEqual( store.lookupInchiKey("UHOVQNZJYSORNB-UHFFFAOYSA-N"), [0])
            self.assertEqual( store.descriptors().get(0), (0.062268570782092456, 0.0, 72.06599999999999, 2.0, 0.062268570782092456, 2.0, -0.062268570782092456, 78.046950192, 78.11399999999999, 30.0, 2.0, 2.0, -0.062268570782092456, 3.000000000000001, 71.96100505779535, 4.242640687119286, 3.464101615137755, 3.464101615137755, 3.0, 2.0000000000000004, 2.0000000000000004, 1.1547005383792521, 1.1547005383792521, 0.6666666666666671, 0.6666666666666671, 0.38490017945975075, 0.38490017945975075, -0.78, 34.3994618804395, 3.4115708812260532, 1.6057694396735218, 0.5823992601400448, 37.43140311949697, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 36.39820241076966, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 36.39820241076966, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 36.39820241076966, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 36.39820241076966, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.6866, 26.441999999999993, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)
                
    def testMany(self):
        try:
            fname = tempfile.mktemp()
            storefname = tempfile.mktemp()
            with open(fname, 'w') as f:
                f.write(many_smiles)
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="RDKit2D",
                                                  index_inchikey=True )
            make_store.make_store(opts)

            store = DescriptaStore(storefname)

            for i in range(10):
                self.assertEqual( store.lookupName(str(i)), i)
                
            for i in range(10):
                m = store.molIndex().getRDMol(i)
                inchi = AllChem.InchiToInchiKey(AllChem.MolToInchi(m))
                self.assertEqual( store.lookupInchiKey(inchi), [i])
                
            self.assertEqual( store.descriptors().get(0), (0.062268570782092456, 0.0, 72.06599999999999, 2.0, 0.062268570782092456, 2.0, -0.062268570782092456, 78.046950192, 78.11399999999999, 30.0, 2.0, 2.0, -0.062268570782092456, 3.000000000000001, 71.96100505779535, 4.242640687119286, 3.464101615137755, 3.464101615137755, 3.0, 2.0000000000000004, 2.0000000000000004, 1.1547005383792521, 1.1547005383792521, 0.6666666666666671, 0.6666666666666671, 0.38490017945975075, 0.38490017945975075, -0.78, 34.3994618804395, 3.4115708812260532, 1.6057694396735218, 0.5823992601400448, 37.43140311949697, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 36.39820241076966, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 36.39820241076966, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 36.39820241076966, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 36.39820241076966, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.6866, 26.441999999999993, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)
