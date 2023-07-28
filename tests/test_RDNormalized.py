from __future__ import print_function
import unittest
from rdkit.Chem import AllChem
import numpy, math
from descriptastorus import make_store, DescriptaStore
from descriptastorus.descriptors import (
    rdDescriptors, rdNormalizedDescriptors, DescriptorGenerator, dists)
import contextlib, tempfile, os, shutil, sys
import datahook

one_smiles = "c1ccccc1 0"
many_smiles = "\n".join( [ "C"*i + "c1ccccc1 " + str(i) for i in range(10) ] )
from expected_normalized_results import expected

def compare_results(unit, result, expected, columns):
    unit.assertEqual(result[0], expected[0])
    for i, (x,y) in enumerate(zip(result[1:], expected[1:])):
        unit.assertAlmostEqual(x,y, 5, msg="At %s"%columns[i+1][0])

class TestCase(unittest.TestCase):
    def testHaveNormalizations(self):
        missing = []
        for feature in rdDescriptors.RDKIT_PROPS[rdDescriptors.CURRENT_VERSION]:
            if feature not in dists.dists:
                missing.append(feature)

        print(repr(missing))
        self.assertEqual(missing, [])

    def testNormalized(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
            with open(fname, 'w') as f:
                f.write(many_smiles),
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="RDKit2DNormalized",
                                                  index_inchikey=True )
            make_store.make_store(opts)
            generator = DescriptorGenerator.REGISTRY["RDKit2DNormalized".lower()]
            results = []
            with contextlib.closing(DescriptaStore(storefname)) as store:
                for i in range(10):
                    r = store.descriptors().get(i)
                    compare_results(self, r, expected[i], generator.GetColumns())

        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)
                
if __name__ == '__main__':  #pragma: no cover
    unittest.main()

