#  Copyright (c) 2018, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
from __future__ import print_function
from . import raw
from . import MolFileIndex
from . import keyvalue
import os, sys, contextlib, pickle
from .keyvalue import KeyValueAPI

import logging
logger = logging.getLogger("descriptastorus")

try:
    from .descriptors import MakeGenerator
except:
    MakeGenerator = None
    logger.exception("Unable to make new descriptors, descriptor generator requirements not properly installed")

from .raw import Mode

class DescriptaStoreIter:
    def __init__(self, store):
        self.store = store
        self.i = -1
    def __next__(self):
        self.i += 1
        if self.i == len(self.store):
            raise StopIteration()

        try:
            return self.store.index.get(self.i), self.store.getDescriptors(self.i)
        except:
            print("== DescriptaStoreIter Failed at index", self.i)
            raise
    def next(self):
        return self.__next__(self)

def get_options(dbdir):
    optionsfile = os.path.join(dbdir, "__options__")
    if os.path.exists(optionsfile):
        with open(optionsfile, 'rb') as f:
            return pickle.load(f)
    else:
        raise IOError("Not a valid Descriptastore, no options file")
    
class DescriptaStore:
    def __init__(self, dbdir, mode=Mode.READONLY):
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
        self.db = raw.RawStore(dbdir, mode=mode)
        self.index = MolFileIndex.MolFileIndex(os.path.join(dbdir, "__molindex__"))
        options = self.options = get_options(dbdir)
        keystore = options.get("keystore", "kyotostore")
        
        key_store_type = None
        if keystore:
            key_store_type = KeyValueAPI.get_store(keystore)
            if not key_store_type:
                logger.warning("Keystore %r not available, skipping", keystore)

        self.inchikey = self.name = None
        if key_store_type:
            inchi = os.path.join(dbdir, key_store_type().get_actual_filename("inchikey"))
            if os.path.exists(inchi):
                self.inchikey = key_store_type()
                self.inchikey.open(os.path.join(dbdir, "inchikey"), mode=mode)

            name = os.path.join(dbdir, key_store_type().get_actual_filename("name"))
            if os.path.exists(name):
                self.name = key_store_type()
                self.name.open(os.path.join(dbdir, "name"), mode=mode)

        # index the calculated flags
        datacols = [(i,name) for i,name in enumerate(self.db.colnames) if "_calculated" not in name]
        self.datanames = [name for i,name in datacols]
        self.dataindices = [i for i,name in datacols]
        

    def close(self):
        self.db.close()
        self.index.close()
        if self.inchikey is not None:
            self.inchikey.close()
        
        if self.name is not None and hasattr(self.name, "close"):
            self.name.close()
            
    def __len__(self):
        return self.db.N

    def __iter__(self):
        return DescriptaStoreIter(self)

    def getDescriptorCalculator(self):
        """Returns the descriptor calculator (if possible) for the store
        In general this requires the same run-time environment as the 
        storage, so this might not be possible"""
        try:
            return MakeGenerator(self.options['descriptors'].split(","))
        except:
            logger.exception("Unable to make generator from store")
            return None

    def getDescriptorNames(self, keepCalculatedFlags=False):
        """keepCalculatedFlags=False -> return the descriptor names for the store
        if keepCalculatedFlags is True return the boolean flags that indicate
        if results were calculated for the descriptor subset.
        """
        if keepCalculatedFlags:
            return self.db.colnames[:]
        return self.datanames
        
    def getDescriptors(self, index, keepCalculatedFlags=False):
        """index, keepCalculatedFlags=False -> return the descriptors at index
        if keepCalculatedFlags is True return the boolean flags that indicate
        if results were calculated for the descriptor subset.
        """
        
        v = self.db.get(index)
        if keepCalculatedFlags:
            return v
        else:
            return [v[i] for i in self.dataindices]
        
    def getDescriptorsAsDict(self, index):
        """index -> return the descriptors as an index"""
        return self.db.getDict(index)
    
    def descriptors(self):
        """Returns the raw storage for the descriptors"""
        return self.db

    def molIndex(self):
        """Returns the mol index"""
        return self.index

    def lookupName(self, name):
        """name -> returns the index of the given name"""
        if self.name is None:
            try:
                logger.warning("Using slower memory intensive option")
                logger.warning("Loading names...")
                self.name = {name:i for i, (moldata, name) in enumerate(self.index)}
                logger.warning("...done loading")
                print(self.name)
            except:
                logger.exception("Names not available from original input")
                raise ValueError("Name index not available")
            assert self.name

        try:
            row = self.name.get(name)
            if row is None:
                raise KeyError(name)
        except:
            raise KeyError("Name %r not found"%name)
        
        return row
    
    def lookupInchiKey(self, key):
        """key -> returns the indicies of the inchi key"""
        if self.inchikey is None:
            raise ValueError("Inchi index not available")
        res = self.inchikey.get(key)
        if res is None:
            raise KeyError(key)
        return res
    
