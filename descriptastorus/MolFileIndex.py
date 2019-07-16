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
from rdkit.Chem import AllChem

import os, numpy, sys
from . import raw
import logging, shutil, pickle
import csv
from io import StringIO
import functools

def SDFNameGetter(buffer):
    return buffer.split("\n")[0].strip()


namefxns = {None: None, "molfile": SDFNameGetter}

def nameOptFile(indexdir):
    return os.path.join(indexdir, "__opts__")

class MolFileIter:
    def __init__(self, raw):
        self.raw = raw
        self.i = -1
        
    def __next__(self):
        self.i += 1
        try:
            if self.i >= len(self.raw):
                raise StopIteration()
            return self.raw.get(self.i)
        except IndexError:
            raise StopIteration()

    def next(self): return self.__next__()

def reader(s, dialect=None):
    return csv.reader(StringIO(s), dialect=dialect)

def whitespace_reader(s):
    return [s.split()]

class MolFileIndex:
    """Index for a molecule file to provide random access to the internal molecules.
    """
    
    def __init__(self, indexDirectory, mode=raw.Mode.READONLY):
        """Fast random access to a smiles file by row index
        indexDirectory = directory of the molfile index

        Example
        -----------
        See MakeSmilesIndex and MakeSDF Index to make indexed files
        """
        optfile = nameOptFile(indexDirectory)

        if os.path.exists(optfile):
            with open(optfile, 'rb') as f:
                options = pickle.load(f)
        else:
            raise IOError("Not a molfile index")


        self.db = raw.RawStore(indexDirectory, mode=mode)
        self._nameGetter = None
        self.filename = options['filename']
        self.hasHeader = options['hasHeader']
        self.smilesColumn = options['smilesColumn']
        self.nameColumn = options['nameColumn']
        self.sep = options['sep']
        self.nameFxnName = options['nameFxnName']

        try:
            if self.sep == None:
                self.sep = ' '
                
            if self.sep == 'excel':
                self.reader = reader
            elif self.sep == 'excel_tab':
                self.reader = functools.partial(reader, dialect=csv.excel_tab)
            elif self.sep == 'whitespace':
                self.reader = whitespace_reader
            else:
                # assume a seperator
                csv.register_dialect('custom_dialect', delimiter=self.sep, skipinitialspace=True)
                self.reader = functools.partial(reader, dialect='custom_dialect')
        except:
            logging.exception("Can't initialize delimiter: %s", self.sep)
            

        self.nameFxn = namefxns[self.nameFxnName]
        
        if self.hasHeader:
            self.N = self.db.N - 3
        else:
            self.N = self.db.N - 2

        # mmap?
        self.filename = os.path.join(indexDirectory, self.filename)        
        self.f = open(self.filename, 'r')
        
        if self.hasHeader:
            colnames = self.colnames = self._get(None)
        else:
            colnames = self.colnames = ["column_%d"%x for x in range(len(self._get(0)))]

        # get the first entry
        row = self._get(0) # takes header into account internally
            
            
        if self.smilesColumn != -1:
            try:
                self.smilesColIdx = int(self.smilesColumn)
            except ValueError:
                try:
                    self.smilesColIdx = colnames.index(self.smilesColumn)
                except ValueError:
                    raise IndexError("Specified smiles column %r name not in header\n"
                                     "\tHeader is %r\n"
                                     "\tPerhaps the seperator is misspecified (currently %r)"%(
                                         self.smilesColumn,
                                         self.colnames,
                                         self.sep)
                                     )
            
            if len(row) <= self.smilesColIdx:
                raise IndexError("Smiles Column %d greater than rowsize %s\n"
                                 "Perhaps the `seperator is mispecified (currently %r)"% (
                                     self.smilesColIdx,
                                     len(row),
                                     self.sep))
        self.nameidx = -1
        if self.nameColumn is not None:
            try:
                self.nameidx = int(self.nameColumn)
            except ValueError:
                try:
                    self.nameidx = colnames.index(self.nameColumn)
                except ValueError:
                    raise IndexError("Specified name column %r name not in header\n"
                                     "\tHeader is %r\n"
                                     "\tPerhaps the seperator is misspecified (currently %r)"%(
                                         self.smilesColumn,
                                         self.colnames,
                                         self.sep)
                                     )


            if len(row) <= self.nameidx:
                raise IndexError("Name Column %d greater than rowsize %s\n"
                                 "Perhaps the seperator is mispecified (currently %r)"% (
                                     self.smilesColIdx,
                                     len(row),
                                     self.sep))

    def __del__(self):
        self.close()
        
    def close(self):
        self.db.close()
        self.f.close()

    def __len__(self):
        return self.N
    
    def __iter__(self):
        return MolFileIter(self)
        
    def _get(self, idx):
        if idx is None:
            idx = 0
        elif self.hasHeader:
            idx += 1
            
        start = self.db.get(idx)[0]
        end = self.db.get(idx+1)[0]
        self.f.seek(start,0)
        buf = self.f.read(end-start-1)
        try:
            if self.smilesColumn != -1:
                return list(self.reader(buf))[0]#buf.split(self.sep)
        except:
            logging.exception("Whups, can't split")
            raise
        return buf

    def header(self):
        """Return header column (throws ValueError if no header column is available)"""
        if self.hasHeader:
            return self._get(None)
        raise ValueError("Datastore doesn't have a header")
    
    def get(self, idx):
        """idx -> gets the data at row idx
        return a list if the data is a smiles like file
        returns a string buffer otherwise
        """
        v = self._get(idx)
        if self.smilesColIdx != -1:
            moldata = v[self.smilesColIdx]
            if self.nameidx != -1:
                name = v[self.nameidx]
                return moldata, name
            return moldata
        if self._nameGetter:
            return v, self._nameGetter(v)
        return v
    
    def getMol(self, idx):
        """Returns input data for the molecule"""
        if self.smilesColIdx != -1:
            return self._get(idx)[self.smilesColIdx]
        return self._get(idx)

    def getRDMol(self, idx):
        """Returns the RDKit molecular representation of the input data"""
        data = self._get(idx)
        if self.smilesColIdx != -1:
            m = AllChem.MolFromSmiles(data[0])
            if m:
                if len(data) > 0:
                    m.SetProp("_Name", data[1])
            return m
        return AllChem.MolFromMolBlock(data)
    
    def getName(self, idx):
        if self.nameidx == -1:
            if self._nameGetter:
                return self._nameGetter(self._get(idx))
            
            raise ValueError("SmilesIndex does not have a name column or a name retriever")
        
        return self._get(idx)[self.nameidx]

def simplecount(filename):
    lines = 0
    with open(filename) as f:
        for line in f:
            lines += 1
    return lines

def index(fname, word):
    fsize = os.path.getsize(fname)
    bsize = 2**16
    
    with open(fname, 'rb') as f:

        buffer = None
        overlap = len(word) - 1
        while True:
            if (f.tell() >= overlap and f.tell() < fsize):
                f.seek(f.tell() - overlap)
            buffer = f.read(bsize)
            if buffer:
                pos = buffer.find(word)
                while pos != -1:
                    yield f.tell() - (len(buffer) - pos)
                    pos = buffer.find(word, pos+1)
                
            else:
                break
    
def MakeSmilesIndex(filename, dbdir, hasHeader, smilesColumn, nameColumn=-1, sep=None,
                    reIndex=False):
    """Make smiles index -> index a smiles file for random access
    filename: filename to index
    dbdir: name of the index store
    hasHeader: do we have a header
    smilesColumn: column for the smiles string
    nameColumn: column for the name string
    sep: seperator for the file i.e. '\t'
    reIndex: reIndex the existing file
    otherwise Copies file over to index"""
    targetFilename = os.path.join(dbdir, os.path.basename(filename))
    if reIndex and os.path.abspath(filename) != os.path.abspath(targetFilename):
        raise ValueError("If reindex is set, filename must be the storage filename (this is a sanity check\n%s\n%s"%(filename, targetFilename))
    
    sz = os.path.getsize(filename)
    
    N = simplecount(filename)

    if sz < 2**8:
        dtype = numpy.uint8
    elif sz < 2**16:
        dtype = numpy.uint16
    elif sz < 2**32:
        dtype = numpy.uint32
    else:
        dtype = numpy.uint64

    db = raw.MakeStore([("index", dtype)], N+2, dbdir,
                       checkDirectoryExists=(not reIndex))
    cpfile = targetFilename
    if not reIndex:
        logging.info("Copying molecule file to index...")
        shutil.copy(filename, cpfile)
        logging.info("Done copying")
    else:
        logging.info("Reindexing existing smiles file...")
        
    options = {'filename': os.path.basename(filename),
               'hasHeader': hasHeader,
               'smilesColumn': smilesColumn,
               'nameColumn': nameColumn,
               'nameFxnName': None,
               'sep': sep}

    # save the options
    optfile = nameOptFile(dbdir)

    with open(optfile, 'wb') as f:
        pickle.dump(options, f)
    
    # first row
    #  TODO sniff newline...
    logging.info("Indexing...")
    db.putRow(0, [0])
    for i,pos in enumerate(index(cpfile, b"\n")):
        db.putRow(i+1, [pos+1])
    db.close()
    return MolFileIndex(dbdir)
#, os.path.basename(filename), smilesColumn,
#                        nameColumn=nameColumn, hasHeader=hasHeader,
#                        sep=sep)

def MakeSDFIndex(filename, dbdir):
    """Make smiles index -> index a smiles file for random access"""
    sadf
    sz = os.path.getsize(filename)
    
    N = simplecount(filename)
    if sz < 2**8:
        dtype = numpy.uint8
    elif sz < 2**16:
        dtype = numpy.uint16
    elif sz < 2**32:
        dtype = numpy.uint32
    else:
        dtype = numpy.uint64

    # TODO sniff newline ...
    indices = list(index(filename, b"$$$$\n"))
    
    db = raw.MakeStore([("index", dtype)], N+1, dbdir)

    # first row
    db.putRow(0, [0])
    for i, idx in enumerate(indices):
        db.putRow(i+1, [pos+1])
    
    return MolFileIndex(filename, dbdir, nameFxn=SDFNameGetter)

        
