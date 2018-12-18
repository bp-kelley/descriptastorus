from ..keyvalue import KeyValueAPI
from ..raw import Mode
import logging, os
import dbm

class DBMStore(KeyValueAPI):
    DB = dbm
    def open(self, fn, mode):
        fn = self.get_actual_filename(fn)

        if mode in [Mode.READONLY, Mode.READONCE]:
            db = self.DB.open(fn, 'r')
        elif mode ==  Mode.WRITE:
            db = self.DB.open(fn, 'c')
        elif mode ==  Mode.APPEND:
            db = self.DB.open(fn, 'w')
        else:
            raise ValueError("Invalide mode %r for opening %s"%(mode, self.__class__.__name__))
        self.db = db
        
    def get_actual_filename(self, fn):
        return fn + ".dumbdb"

    def close(self):
        self.db.close()

    def get_raw(self, key):
        return self.db[key]

    def set_raw(self, key, value):
        self.db[key] = value

    def __contains__(self, key):
        return key in self.db


KeyValueAPI.register("dbmstore", DBMStore)


