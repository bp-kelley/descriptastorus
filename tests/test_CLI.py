from __future__ import print_function
import unittest
from descriptastorus import append_store, make_store
from descriptastorus.cli import storus
from descriptastorus.DescriptaStore import DescriptaStore
from descriptastorus.descriptors.rdDescriptors import RDKit2D
import tempfile, contextlib
import logging, os, shutil
from rdkit.Chem import AllChem

logging.getLogger().setLevel(logging.INFO)
make_store.DEFAULT_KEYSTORE = "dbmstore"
one_smiles = "c1ccccc1 0"
many_smiles = "\n".join( [ "C"*i + "c1ccccc1 " + str(i) for i in range(10) ] + ["NOSTRUCT foo"] )
many_smiles2 = "\n".join( [ "C"*i + "c1ccccc1 " + str(i+11) for i in range(10) ] + ["NOSTRUCT foo2"] )

class RDKit2DSubset(RDKit2D):
    NAME="RDKit2DSubset"
    def __init__(self):
        RDKit2D.__init__(self, properties=[
            'ExactMolWt',
            'NumAliphaticRings', 'NumAromaticCarbocycles',
            'NumAromaticHeterocycles', 'NumAromaticRings'])
RDKit2DSubset()

from descriptastorus.descriptors import MakeGenerator
MakeGenerator(['rdkit2dsubset'])

def toDict( v ):
    return {n:x for n,x in zip([
            'RDKit2DSubset_calculated', 'ExactMolWt',
            'NumAliphaticRings', 'NumAromaticCarbocycles',
            'NumAromaticHeterocycles', 'NumAromaticRings'],
                               v)}

class TestCase(unittest.TestCase):
    def setUp(self):
        RDKit2DSubset()
        
    def testMakeStore(self):
        print("*"*44)
        fname = tempfile.mktemp()+".smi"
        storefname = tempfile.mktemp()+".store"
        with open(fname, 'w') as f:
            f.write(many_smiles)

        args = ["--index-inchikey",
                "--smilesColumn",  "0",
                "--nameColumn", "1",
                "--seperator",  " ",
                "--numprocs",  "1",
                "--descriptors", "RDKit2DSubset",
                fname, storefname]

        try:
            opts = storus.parser.parse_args(args)
            logging.error(repr(opts))
            make_store.make_store(make_store.MakeStorageOptions(**vars(opts)))
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

                

    def testAppendStore(self):
        fname = tempfile.mktemp()+".smi"
        fname2 = tempfile.mktemp()+"2.smi"
        storefname = tempfile.mktemp()+".store"
        with open(fname, 'w') as f:
            f.write(many_smiles)
        with open(fname2, 'w') as f:
            f.write(many_smiles2)

        # make the first store
        args = ["--index-inchikey",
                "--smilesColumn",  "0",
                "--nameColumn", "1",
                "--seperator",  " ",
                "--numprocs",  "1",
                "--descriptors", "RDKit2DSubset",
                fname, storefname]

        args2 = ["--append",
                 "--index-inchikey",
                 "--smilesColumn",  "0",
                 "--nameColumn", "1",
                 "--seperator",  " ",
                 "--numprocs",  "1",
                 "--descriptors", "RDKit2DSubset",
                 fname2, storefname]

        try:
            opts = storus.parser.parse_args(args)
            make_store.make_store(make_store.MakeStorageOptions(**vars(opts)))

            opts = storus.parser.parse_args(args2)
            append_store.append_smiles(append_store.AppendStorageOptions(**vars(opts)))

            with contextlib.closing(DescriptaStore(storefname)) as store:
                assert len(store) == 22, str(len(store))
                for i in range(20):
                    m = store.molIndex().getRDMol(i)
                    if m:
                        sm = AllChem.MolToSmiles(m)
                        inchi = AllChem.InchiToInchiKey(AllChem.MolToInchi(m))
                        logging.info("%s: %s %s"%(i, sm, inchi))
                    else:
                        logging.info("%s: nostruct"%i)
                    

                for i in range(10):
                    self.assertEqual( store.lookupName(str(i)), i)
                for i in range(10):
                    m = store.molIndex().getRDMol(i)
                    inchi = AllChem.InchiToInchiKey(AllChem.MolToInchi(m))
                    self.assertEqual( store.lookupInchiKey(inchi), [i, i+11])

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
                    if m:
                        sm = AllChem.MolToSmiles(m)
                        inchi = AllChem.InchiToInchiKey(AllChem.MolToInchi(m))
                        logging.info("%s: %s"%(i, inchi))
                    
                        self.assertEqual( store.lookupInchiKey(inchi), [i, i+11])
                    v = store.descriptors().get(i)
                    sv = tuple(calc.process(sm))
                    self.assertEqual(v, sv)

                    
                for i in range(10):
                    m = store.molIndex().getRDMol(i)
                    if not m: continue
                    
                    inchi = AllChem.InchiToInchiKey(AllChem.MolToInchi(m))
                    m = store.molIndex().getRDMol(i+11)
                    self.assertTrue(m!=None)
                    inchi2 = AllChem.InchiToInchiKey(AllChem.MolToInchi(m))
                    self.assertEqual(inchi, inchi2)
                    self.assertEqual( store.lookupInchiKey(inchi), [i, i+11])
                    
                for i in range(2):
                    self.assertEqual(store.descriptors().get(11+0), (True, 78.046950192, 0.0, 1.0, 0.0, 1.0))
                    self.assertEqual(store.descriptors().get(11+1), (True, 92.062600256, 0.0, 1.0, 0.0, 1.0))
                    self.assertEqual(store.descriptors().get(11+2), (True, 106.07825032, 0.0, 1.0, 0.0, 1.0))
                    self.assertEqual(store.descriptors().get(11+3), (True, 120.093900384, 0.0, 1.0, 0.0, 1.0))
                    self.assertEqual(store.descriptors().get(11+4), (True, 134.109550448, 0.0, 1.0, 0.0, 1.0))
                    self.assertEqual(store.descriptors().get(11+5), (True, 148.125200512, 0.0, 1.0, 0.0, 1.0))
                    self.assertEqual(store.descriptors().get(11+6), (True, 162.140850576, 0.0, 1.0, 0.0, 1.0))
                    self.assertEqual(store.descriptors().get(11+7), (True, 176.15650064, 0.0, 1.0, 0.0, 1.0))
                    self.assertEqual(store.descriptors().get(11+8), (True, 190.172150704, 0.0, 1.0, 0.0, 1.0))
                    self.assertEqual(store.descriptors().get(11+9), (True, 204.187800768, 0.0, 1.0, 0.0, 1.0))
                    self.assertEqual(store.descriptors().get(11+10), (False, 0.0, 0.0, 0.0, 0.0, 0.0))                       
        finally:
            if os.path.exists(fname):
                os.unlink(fname)
            if os.path.exists(storefname):
                shutil.rmtree(storefname)
if __name__ == '__main__':  #pragma: no cover
    unittest.main()
