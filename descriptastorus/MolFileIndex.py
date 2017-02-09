from __future__ import print_function
from rdkit.Chem import AllChem

import os, numpy, sys
from . import raw
import shutil, pickle

def SDFNameGetter(buffer):
    return buffer.split("\n")[0].strip()


namefxns = {None: None, "molfile": SDFNameGetter}

def nameOptFile(indexdir):
    return os.path.join(indexdir, "__opts__")

class MolFileIter:
    def __init__(self, raw):
        self.raw = raw
        self.i = -1
    def next(self):
        self.i += 1
        if self.i >= len(self.raw):
            raise StopIteration()
        return self.raw.get(self.i)

class MolFileIndex:
    """Index for a molecule file to provide random access to the internal molecules.
    """
    
    def __init__(self, indexDirectory):
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


        self.db = raw.RawStore(indexDirectory)
        self._nameGetter = None
        self.filename = options['filename']
        self.hasHeader = options['hasHeader']
        self.smilesColumn = options['smilesColumn']
        self.nameColumn = options['nameColumn']
        self.sep = options['sep']
        self.nameFxnName = options['nameFxnName']
        
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
                                 "Perhaps the seperator is mispecified (currently %r)"% (
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
    def close(self):
        self.db.close()
        self.f.close()

    def __len__(self):
        return self.db.N
    
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
        if self.smilesColumn != -1:
            return buf.split(self.sep)
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
    
def MakeSmilesIndex(filename, dbdir, hasHeader, smilesColumn, nameColumn=-1, sep=None):
    """Make smiles index -> index a smiles file for random access
    n.b. Copies file over to index"""
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

    db = raw.MakeStore([("index", dtype)], N+2, dbdir)
    cpfile = os.path.join(dbdir, os.path.basename(filename))
    print("Copying molecule file to index...", file=sys.stderr)
    shutil.copy(filename, cpfile)
    print("Done copying", file=sys.stderr)
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
    print("Indexing...", file=sys.stderr)
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

        
