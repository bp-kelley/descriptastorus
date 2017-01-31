from __future__ import print_function
from . import raw
from . import MolFileIndex
import os, sys, contextlib, pickle

try:
    import kyotocabinet
except ImportError:
    kyotocabinet = None

class DescriptaStoreIter:
    def __init__(self, store):
        self.store = store
        self.i = -1
    def next(self):
        self.i += 1
        if self.i == len(self.store):
            raise StopIteration()

        try:
            return self.store.index.get(self.i), self.store.db.get(self.i)
        except:
            print("== DescriptaStoreIter Failed at index", self.i)
            raise
    
class DescriptaStore:
    def __init__(self, dbdir):
        """dbdir -> opens a descriptor storage
         
        >>> store = DescriptaStore(db)
        >>> len(store)

        # access the options used to create this store
        #  (this is optional and may not exist)
        >>> store.options
        ...
        
        Iterate through molecule data ([moldata, <optional name>], descriptors)
        >>> for moldata, descriptors in store:
        >>>     pass

        Iterate through only the descriptors
        >>> for i,prop in enumerate(store.descriptors()):
        >>>    pass

        If name indexed:
        >>> row = store.lookupName("ZWIMER-03065")
        
        If inchi key index:
        Since inchi keys may collide, this can return multiple indices
        >>>  rows = store.lookupInchiKey("BCWYEXBNOWJQJV-UHFFFAOYSA-N")
        """
        self.desctiporDB = dbdir
        self.db = raw.RawStore(dbdir)
        self.index = MolFileIndex.MolFileIndex(os.path.join(dbdir, "__molindex__"))

        inchi = os.path.join(dbdir, "inchikey.kch")
        if os.path.exists(inchi):
            if not kyotocabinet:
                print("Inchi lookup exists, but kyotocabinet is not installed.",
                      file=sys.stderr)
            else:
                self.inchikey = kyotocabinet.DB()
                self.inchikey.open(inchi, kyotocabinet.DB.OREADER)
        else:
            self.inchikey = None

        name = os.path.join(dbdir, "name.kch")
        if os.path.exists(name):
            if not kyotocabinet:
                print("Name lookup exists, but kyotocabinet is not installed.",
                      file=sys.stderr)
            else:
                self.name = kyotocabinet.DB()
                self.name.open(name, kyotocabinet.DB.OREADER)
        else:
            print("Couldn't open name db", name, file=sys.stderr)
            self.name = None

        self.options = None
        optionsfile = os.path.join(dbdir, "__options__")
        if os.path.exists(optionsfile):
            with open(optionsfile) as f:
                self.options = pickle.load(f)

    def close(self):
        self.db.close()
        self.index.close()
        if self.inchikey is not None:
            self.inchikey.close()
        
        if self.name is not None:
            self.name.close()
            
    def __len__(self):
        return self.db.N

    def __iter__(self):
        return DescriptaStoreIter(self)
    
    def descriptors(self):
        """Returns the raw storage for the descriptors"""
        return self.db

    def molIndex(self):
        """Returns the mol index"""
        return self.index

    def lookupName(self, name):
        """name -> returns the index of the given name"""
        if self.name is None:
            raise ValueError("Name index not available")

        try:
            row = int(self.name[name])
        except:
            raise IndexError("Name %r not found"%name)
        
        return row
    
    def lookupInchiKey(self, key):
        """key -> returns the indicies of the inchi key"""
        if self.inchikey is None:
            raise ValueError("Inchi index not available")
        res =  eval(self.inchikey[key])
        return res
    
