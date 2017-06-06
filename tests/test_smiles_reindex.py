from __future__ import print_function

import unittest
from descriptastorus import MolFileIndex
from descriptastorus import append_store
import os, shutil, sys
import logging

import datahook


TEST_DIR = "test-smiles-reindex"

class TestCase(unittest.TestCase):
    def setUp(self):
        if os.path.exists(TEST_DIR):
            shutil.rmtree(TEST_DIR, ignore_errors=True)
    def tearDown(self):
        if os.path.exists(TEST_DIR):
            shutil.rmtree(TEST_DIR, ignore_errors=True)

            
    def testReIndexing(self):
        index = MolFileIndex.MakeSmilesIndex(
            os.path.join(datahook.datadir, "../data/test1.smi"), TEST_DIR, hasHeader=True,
            smilesColumn="smiles", nameColumn="name")
        
        
        logging.info("Running index test")
        self.assertEqual(index.N, 14)
        self.assertEqual(index.getMol(12), 'c1ccccc1CCCCCCCCCCCC')
        self.assertEqual(index.getName(12), '13')

        self.assertEqual(index.getRDMol(13), None)

        index.close()
        fname = index.filename
        # just append the same file
        append_store.append_smiles_file(os.path.join(datahook.datadir, "../data/test1.smi"),
                                        fname, hasHeader=True)
        # now reindex
        logging.debug(fname)
        with open(fname) as f:
            logging.debug(f.read())
            
        index =  MolFileIndex.MakeSmilesIndex(
            fname, TEST_DIR, hasHeader=True,
            smilesColumn="smiles", nameColumn="name", reIndex=True)

        self.assertEqual(index.N, 28)
        self.assertEqual(index.getMol(12), 'c1ccccc1CCCCCCCCCCCC')
        self.assertEqual(index.getName(12), '13')

        self.assertEqual(index.getMol(12+14), 'c1ccccc1CCCCCCCCCCCC')
        self.assertEqual(index.getName(12+14), '13')
        index.close()
        
    def testReIndexingNoLineEnding(self):
        index = MolFileIndex.MakeSmilesIndex(
            os.path.join(datahook.datadir, "../data/test2-nolineending.smi"),
            TEST_DIR, hasHeader=True,
            smilesColumn="smiles", nameColumn="name")
        
        
        logging.info("Running index test")
        self.assertEqual(index.N, 14)
        self.assertEqual(index.getMol(12), 'c1ccccc1CCCCCCCCCCCC')
        self.assertEqual(index.getName(12), '13')

        self.assertEqual(index.getRDMol(13), None)

        index.close()
        fname = index.filename
        # just append the same file
        append_store.append_smiles_file(os.path.join(datahook.datadir, "../data/test1.smi"),
                                        fname, hasHeader=True)
        # now reindex
        index =  MolFileIndex.MakeSmilesIndex(
            fname, TEST_DIR, hasHeader=True,
            smilesColumn="smiles", nameColumn="name", reIndex=True)

        self.assertEqual(index.N, 28)
        self.assertEqual(index.getMol(12), 'c1ccccc1CCCCCCCCCCCC')
        self.assertEqual(index.getName(12), '13')

        self.assertEqual(index.getMol(12+14), 'c1ccccc1CCCCCCCCCCCC')
        self.assertEqual(index.getName(12+14), '13')


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    unittest.main()
