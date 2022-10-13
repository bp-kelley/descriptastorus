from __future__ import print_function
import sys, unittest
from descriptastorus import keyvalue
from descriptastorus.mode import Mode
import logging, os, shutil, numpy, random, time

try:
    import kyotocabinet
except:
    kyotocabinet = None

class TestKeyValue(unittest.TestCase):
    directory = "mydir5"

    def setUp(self):
        if os.path.exists(self.directory):
            logging.info("startup: Removing directory")
            shutil.rmtree(self.directory)
        os.mkdir(self.directory)
        
    def tearDown(self):
        if os.path.exists(self.directory):
            logging.info("tear down: Removing directory")
            shutil.rmtree(self.directory)

    @unittest.skipUnless(kyotocabinet, "kyotocabinet is not available")
    def test_set_get_kyoto(self):
        store = keyvalue.KeyValueAPI.get_store("kyotostore")
        if store:
            s = store()
            self.assertTrue(os.path.exists(self.directory))
            fn = os.path.join(self.directory, "foo")
            s.open(fn, Mode.WRITE)
            s.set('0', 1)
            self.assertEqual(s.get('0'), 1)

            s.set('1', [1,2,3,4])
            self.assertEqual(s.get('1'), [1,2,3,4])
            
    def test_set_get_dbm(self):
        store = keyvalue.KeyValueAPI.get_store("dbmstore")

        s = store()
        self.assertTrue(os.path.exists(self.directory))
        fn = os.path.join(self.directory, "foo")
        s.open(fn, Mode.WRITE)
        s.set('0', 1)
        self.assertEqual(s.get('0'), 1)
        
        s.set('1', [1,2,3,4])
        self.assertEqual(s.get('1'), [1,2,3,4])
        s.close()

        s.open(fn, Mode.READONLY)
        self.assertEqual(s.get('1'), [1,2,3,4])
        s.close()

        

if __name__ == '__main__':  #pragma: no cover
    unittest.main()
        
        
