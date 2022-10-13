from __future__ import print_function
import unittest
from rdkit.Chem import AllChem
from descriptastorus import append_store, make_store, DescriptaStore
from descriptastorus.descriptors.rdDescriptors import RDKit2D
from descriptastorus import descriptors
import logging, struct
import contextlib, tempfile, os, shutil, sys, shutil
import datahook
make_store.DEFAULT_KEYSTORE = "dbmstore"

logging.getLogger().setLevel(logging.DEBUG)

one_smiles = "c1ccccc1 0"
many_smiles = "\n".join( [ "C"*i + "c1ccccc1 " + str(i) for i in range(10) ] + ["NOSTRUCT foo"] )
many_smiles2 = "\n".join( [ "C"*i + "c1ccccc1 " + str(i+11) for i in range(10) ] + ["NOSTRUCT foo2"] )

from rdkit.Chem import Descriptors

class RDKit2DSubset(RDKit2D):
    NAME="RDKit2DSubset"
    def __init__(self):
        RDKit2D.__init__(self, properties=[
            'ExactMolWt',
            'NumAliphaticRings', 'NumAromaticCarbocycles',
            'NumAromaticHeterocycles', 'NumAromaticRings'])
RDKit2DSubset()

class RDKit2DSubsetSmall(RDKit2D):
    NAME="RDKit2DSubsetSmall"
    def __init__(self):
        RDKit2D.__init__(self, properties=[
            'ExactMolWt',])

RDKit2DSubsetSmall()

def toDict( v ):
    return {n:x for n,x in zip([
            'RDKit2DSubset_calculated', 'ExactMolWt',
            'NumAliphaticRings', 'NumAromaticCarbocycles',
            'NumAromaticHeterocycles', 'NumAromaticRings'],
                               v)}

class TestCase(unittest.TestCase):
    def testGenerator(self):
        gen = descriptors.MakeGenerator(("RDKit2DSubset",))
        self.assertTrue(gen!=None)
        
    def testOffByOne(self):
        try:
            fname = tempfile.mktemp()+".smi"
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

    def testColCache(self):
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
                cols = []
                # get normal data
                for idx,_ in enumerate(store.db.colnames):
                    col = list(store.db.getColByIdx(idx))
                    cols.append(col)

                # cache the columns
                store.db.cacheColumns()
                # make sure the datafiles are written
                for idx, _ in enumerate(store.db.colnames):
                    fn=os.path.join(store.db.colCacheDir, str(idx))
                    self.assertTrue(os.path.exists(fn), fn)
                    
                for idx,_ in enumerate(store.db.colnames):
                    col = list(store.db.getColByIdx(idx))
                    self.assertEqual(col,cols[idx])

                # swap a data file
                idx = 0,1
                fn0=os.path.join(store.db.colCacheDir, str(0))
                fn1=os.path.join(store.db.colCacheDir, str(1))
                shutil.move(fn0, fn0+".bak")
                shutil.move(fn1, fn0)
                shutil.move(fn0+".bak", fn1)

                try:
                    col = list(store.db.getColByIdx(0))
                    caught = false
                except struct.error:
                    caught = True

                self.assertTrue(caught, "moving cache file should have broken the cache")

                
                    
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
                
    def testAppend(self):
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

                for i in range(10):
                    m = store.molIndex().getRDMol(i)
                    inchi = AllChem.InchiToInchiKey(AllChem.MolToInchi(m))
                    self.assertEqual( store.lookupInchiKey(inchi), [i])
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

            fname = tempfile.mktemp()+".smi"
            with open(fname, 'w') as f:
                f.write(many_smiles2)
                
            opts.smilesfile = fname
            append_store.append_smiles(opts)
            with contextlib.closing(DescriptaStore(storefname)) as store:
                for i in range(10):
                    self.assertEqual( store.lookupName(str(i)), i)
                    
                for i in range(10):
                    m = store.molIndex().getRDMol(i)
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

    def testAppendStore(self):
        try:
            fname = tempfile.mktemp()+".smi"
            storefname = tempfile.mktemp()+".store"
            storefname2 = tempfile.mktemp()+".store"
            with open(fname, 'w') as f:
                f.write(many_smiles)

            fname2 = tempfile.mktemp()+".smi"
            with open(fname2, 'w') as f:
                f.write(many_smiles2)
                
            opts = make_store.MakeStorageOptions( storage=storefname, smilesfile=fname,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="RDKit2DSubset",
                                                  index_inchikey=True )
            make_store.make_store(opts)

            opts = make_store.MakeStorageOptions( storage=storefname2, smilesfile=fname2,
                                                  hasHeader=False,
                                                  smilesColumn=0, nameColumn=1,
                                                  seperator=" ", descriptors="RDKit2DSubset",
                                                  index_inchikey=True )
            make_store.make_store(opts)
            
            with contextlib.closing(DescriptaStore(storefname)) as store:

                for i in range(10):
                    self.assertEqual( store.lookupName(str(i)), i)

                for i in range(10):
                    m = store.molIndex().getRDMol(i)
                    inchi = AllChem.InchiToInchiKey(AllChem.MolToInchi(m))
                    self.assertEqual( store.lookupInchiKey(inchi), [i])
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

                
            opts.smilesfile = storefname2
            opts.storage = storefname
            append_store.append_store(opts)

            with contextlib.closing(DescriptaStore(storefname)) as store:
                #for i in range(10):
                #    self.assertEqual( store.lookupName(str(i)), i)

                #for i in range(10):
                #    m = store.molIndex().getRDMol(i)
                #    inchi = AllChem.InchiToInchiKey(AllChem.MolToInchi(m))
                #    self.assertEqual( store.lookupInchiKey(inchi), [i])
                    
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
