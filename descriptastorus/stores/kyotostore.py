from ..keyvalue import KeyValueAPI
from ..raw import Mode
import logging, os

logger = logging.getLogger("descriptastorus")

try:
    import kyotocabinet

    class KyotoStore(KeyValueAPI):
        STORE = "kyotocabinet"
        def open(self, fn, mode):
            fn = self.get_actual_filename(fn)

            cabinet = self.cabinet = kyotocabinet.DB()
            if mode in [Mode.READONLY, Mode.READONCE]:
                flags = kyotocabinet.DB.OREADER
            elif mode ==  Mode.WRITE:
                flags = kyotocabinet.DB.OWRITER | kyotocabinet.DB.OCREATE
            elif mode ==  Mode.APPEND:
                flags = kyotocabinet.DB.OWRITER
            else:
                raise ValueError("Invalide mode %r for opening %s"%(mode, self.__class__.__name__))
            cabinet.open(fn, flags)
            
        def close(self):
            self.cabinet.close()
            
        def get_actual_filename(self, fn):
            fn += ".kch"
            return fn
        
        def get_raw(self, key):
            return self.cabinet[key]

        def set_raw(self, key, value):
            self.cabinet[key] = value

        def __contains__(self, key):
            return key in self.cabinet
        
    KeyValueAPI.register("kyotostore", KyotoStore)
except ImportError:
    pass



