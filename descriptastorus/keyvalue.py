from .raw import Mode
import logging
    
class KeyValueAPI:
    """Simple API to wrap various key value stores"""
    REGISTRY = {}
    
    @staticmethod
    def register(name, store):
        KeyValueAPI.REGISTRY[name] = store

    @staticmethod
    def get_store(name):
        res = KeyValueAPI.REGISTRY.get(name, None)
        if res is None:
            logging.warning("Failed to retrieve key value store type %r", name)
        return res

    def get_actual_filename(self, filename):
        raise NotImplementedError
    
    
    def open(self, filename, mode=Mode.READONLY):
        raise NotImplementedError

    def close(self):
        raise NotImplementedError
    
    def get_raw(self, key):
        """Get the value for a particular key (str or bytes
        depending on the python version)"""
        raise NotImplementedError

    def set_raw(self, key, value):
        """Set the value for a particular key.  Values 
        must be str or bytes depending on the python version"""
        raise NotImplementedError

    def get(self, key, default=None):
        try:
            return eval(self.get_raw(key))
        except:
            return default

    def set(self, key, value):
        self.set_raw(key, repr(value))

    
    def __contains__(self, k):
        return NotImplementedError
    
