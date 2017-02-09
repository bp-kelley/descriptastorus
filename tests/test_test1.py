import unittest
from descriptastorus import MolFileIndex
import os, shutil
import logging

import datahook

TEST_DIR = "test1"

class TestCase(unittest.TestCase):
    def setUp(self):
        if os.path.exists(TEST_DIR):
            shutil.rmtree(TEST_DIR, ignore_errors=True)
        index = self.index = MolFileIndex.MakeSmilesIndex(
            os.path.join(datahook.datadir, "../data/test1.smi"), TEST_DIR, hasHeader=True,
            smilesColumn="smiles", nameColumn="name")
        
    def tearDown(self):
        if os.path.exists(TEST_DIR):
            shutil.rmtree(TEST_DIR, ignore_errors=True)

            
    def testIndexing(self):
        logging.info("Running index test")
        self.assertEqual(self.index.N, 14)
        self.assertEqual(self.index.getMol(12), 'c1ccccc1CCCCCCCCCCCC')
        self.assertEqual(self.index.getName(12), '13')

        self.assertEqual(self.index.getRDMol(13), None)


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    unittest.main()
