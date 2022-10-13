from ..keyvalue import KeyValueAPI
from ..raw import Mode
import logging, os
import dbm

class DBMStore(KeyValueAPI):
    DB = dbm
    STORE = "dbm"
    def open(self, fn, mode):
        filename = self._get_dbm_name(fn)

        if mode in [Mode.READONLY, Mode.READONCE]:
            db = self.DB.open(filename, 'r')
        elif mode ==  Mode.WRITE:
            db = self.DB.open(filename, 'c')
        elif mode ==  Mode.APPEND:
            db = self.DB.open(filename, 'w')
        else:
            raise ValueError("Invalid mode %r for opening %s"%(mode, self.__class__.__name__))

        assert os.path.exists(self.get_actual_filename(fn)), self.get_actual_filename(fn)
        self.db = db

    def _get_dbm_name(self, fn):
        return fn + ".dumbdb"
    
    def get_actual_filename(self, fn):
        return self._get_dbm_name(fn) + ".db"

    def close(self):
        self.db.close()

    def get_raw(self, key):
        return self.db[key]

    def set_raw(self, key, value):
        self.db[key] = value

    def __contains__(self, key):
        return key in self.db


KeyValueAPI.register("dbmstore", DBMStore)


