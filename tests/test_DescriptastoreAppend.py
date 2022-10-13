from __future__ import print_function
import unittest
from rdkit.Chem import AllChem
from descriptastorus import append_store, make_store, DescriptaStore
from descriptastorus.descriptors.rdDescriptors import RDKit2D

import contextlib, tempfile, os, shutil, sys
import datahook

make_store.DEFAULT_KEYSTORE = "dbmstore"

one_smiles = "c1ccccc1 0"
two_smiles = "c1ccccc1 1"

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
    def testAppend(self):
        try:
            fname = tempfile.mktemp()+".smi"
            fname2 = tempfile.mktemp()+"-2.smi"
            
            storefname = tempfile.mktemp()+".store"
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

            # now append some junk


            with open(fname2, 'w') as f:
                f.write(two_smiles)

            opts.smilesfile = fname2
            append_store.append_smiles(opts)
            with contextlib.closing(DescriptaStore(storefname)) as store:
                self.assertEqual(len(store), 2)
                self.assertEqual( store.lookupName("0"), 0)

                self.assertEqual( store.lookupInchiKey("UHOVQNZJYSORNB-UHFFFAOYSA-N"), [0,1])
                self.assertEqual(store.descriptors().get(0), (True, 78.046950192, 0.0, 1.0, 0.0, 1.0))
                self.assertEqual( store.lookupName("1"), 1)

                self.assertEqual( store.lookupInchiKey("UHOVQNZJYSORNB-UHFFFAOYSA-N"), [0,1])
                self.assertEqual(store.descriptors().get(1), (True, 78.046950192, 0.0, 1.0, 0.0, 1.0))
            

        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(fname2):
                os.unlink(fname2)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)
                
if __name__ == '__main__':  #pragma: no cover
    unittest.main()
