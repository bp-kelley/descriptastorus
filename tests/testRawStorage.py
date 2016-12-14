import unittest
from descriptastorus import raw
from descriptastorus.descriptors import rdDescriptors
from descriptastorus.descriptors import MakeGenerator

import os, shutil

class TestCase(unittest.TestCase):
    def test_raw(self):
        try:
            directory = "test-store"
            if os.path.exists(directory):
                shutil.rmtree(directory)

            descriptors = MakeGenerator("Morgan3Counts,RDKit2D".split(","))
            N=10
            r = raw.MakeStore(descriptors.GetColumns(), N, directory)

            smiles = ["C"]
            data = []
            for i in range(N):
                row = descriptors.process( "".join(smiles) )
                data.append(row)
                r.putRow( i, row )
                smiles.append("C")

            try:
                r.putRow( i+1, descriptors.process( "".join(smiles) ) )
                self.assertFalse("Shouldn't be able to write past the end of the db")

            except:
                pass

            for i, row in enumerate(data):
                self.assertEquals(tuple(row), r.get(i))
        finally:
            if os.path.exists(directory):
                shutil.rmtree(directory)
                      
if __name__ == '__main__':
    unittest.main()
