"""Raw storage of integer and floating point data.

This provides fast random access to indexed data.

r = RawStore(directory)
row = r.get(0)
row = r.get(100000)

Data can be updated using putRow when opened in append mode.
"""
from __future__ import print_function
import pickle, numpy, os, mmap, struct, sys

# raw stores are little endian!

class Mode:
    READONLY = 0
    WRITE = 1
    APPEND = 2

class RawStoreIter:
    def __init__(self, raw):
        self.raw = raw
        self.i = -1
    def next(self):
        self.i += 1
        if self.i == self.raw.N:
            raise StopIteration()
        return self.raw.get(self.i)
    
class RawStore:
    def __init__(self, directory, mode=Mode.READONLY):
        """Raw storage engine
        directory = existing directory to read or prepare the raw storage
        if N is None, an existing directory is created with N entries"""
        if not os.path.exists(directory):
            raise IOError("Directory %s for raw store does not exist"%
                          directory)
        with open(os.path.join(directory, "__rawformat__"), 'rb') as rawformat:
            self.__dict__.update(pickle.load(rawformat))
            
        fname = self.fname = os.path.join(directory, "__store___")

        if mode == Mode.APPEND:
            self.f = open(fname, 'wb')
        elif mode == Mode.READONLY:
            self._f = open(fname, 'rb')
            access = mmap.ACCESS_READ
        else:
            self._f = open(fname, 'r+b')
            access = mmap.ACCESS_WRITE
            
        self.f = mmap.mmap(self._f.fileno(), 0, access=access)

    def close(self):
        self.f.close()
        self._f.close()

    def __len__(self):
        return self.N
    
    def __iter__(self):
        return RawStoreIter(self)
    
    def get(self, idx):
        """Return the row at idx"""
        offset = idx * self.rowbytes
        try:
            self.f.seek(offset,0)
        except ValueError:
            print("Could not seek to index %d at offset %f offset"%(idx, offset), file=sys.stderr)
            raise IndexError("out or range %d"%idx)
        bytes = self.f.read(self.rowbytes)
        return struct.unpack(self.pack_format, bytes)

    def getColByIdx(self, column):
        """Return the data in the entire column (lazy generator)"""
        # figure out the column
        
        # figure out how many bytes to skip per row
        skip_format = self.pack_format[:column]
        offset = len(struct.pack(skip_format,
                                 *[dtype(0)
                                   for dtype in self.dtypes[:column]]))


        # figure out the bytes to read via the format for the row
        pack_format = self.pack_format[column]
        nbytes = len(struct.pack(pack_format, self.dtypes[column](0)))

        skip = offset - nbytes

        # seek to the first column
        self.f.seek(offset,0)
        # compute the next entry
        offset += self.rowbytes
        while 1:
            bytes = self.f.read(nbytes)
            if bytes == '':
                raise StopIteration
            try:
                yield struct.unpack(pack_format, bytes)
            except struct.error:
                if len(bytes) != nbytes:
                    raise StopIteration
                raise

            try:
                self.f.seek(offset, 0)
            except ValueError:
                raise StopIteration

            
            offset += self.rowbytes
            
    def getCol(self, column_name):
        """Return the column addressed by column_name
        throws IndexError if the column doesn't exist."""
        idx = self.colnames.index(column_name)
        return self.getColByIdx(idx)
        
    def putRow(self, idx, row):
        """Put data row in to row idx.
        Checks to see if the data in v is compatible with the column
        formats.  Throws ValueError on failure."""
        if idx >= self.N:
            raise IndexError("Attempting to write index %s, raw store only has %d rows"%(
                idx, self.N))
        
        if len(row) != len(self.dtypes):
            raise ValueError(
                "data value only has %s entries, "
                "raw store has %s columns"%(
                    len(row), len(self.dtypes)))

        # checks row datatypes
        v = [dtype(v) for dtype,v in zip(self.dtypes, row)]
        offset = idx * self.rowbytes
        self.f.seek(offset,0)
        try:
            bytes = struct.pack(self.pack_format, *row)
        except struct.error:
            print("Can't write row %r"%(row), file=sys.stderr)
            raise
        self.f.write(bytes)

    def write(self, row):
        """Writes row with no datatype checking to the end of the file.
        Do not use with putRow/getRow
        Rows must be written in the correct order"""
        bytes = struct.pack(self.pack_format, *row)
        self.f.write(bytes)

def MakeStore(cols, N, directory, checkDirectoryExists=True):
    if not os.path.exists(directory):
        os.mkdir(directory)
    else:
        if checkDirectoryExists:
            raise IOError("Directory %r for raw store already exists"%directory)

    if not N:
        raise ValueError(
            "When creating a RawStore the total number of elements "
            "must be known,\nplease specify N=?")


    types = []
    names = []
    dtypes = []
    for name, dtype in cols:
        names.append(name)
        if dtype == numpy.int32:
            type = "i"
            dtypes.append(int)
        elif dtype == numpy.int64:
            type = "q"
            dtypes.append(int)
        elif dtype == numpy.uint8:
            type = "B"
            dtypes.append(int)
        elif dtype == numpy.uint16:
            type = "H"
            dtypes.append(int)
        elif dtype == numpy.uint32:
            type = "I"
            dtypes.append(int)
        elif dtype == numpy.uint64:
            type = "Q"
            dtypes.append(int)
        elif dtype == numpy.float32:
            type = "f"
            dtypes.append(float)
        elif dtype == numpy.float64:
            type = "d"
            dtypes.append(float)
        elif dtype == numpy.bool:
            type = "?"
            dtypes.append(bool)
        else:
            raise ValueError("Unhandled numpy type %s"%dtype)

        types.append(type)
                     

    pack_format = "".join(types)
    rowbytes = len(struct.pack(pack_format,
                               *[d(0) for d in dtypes]))
    with open(os.path.join(directory, "__rawformat__"), 'wb') as rawformat:
        pickle.dump(
            {'rowbytes':rowbytes,
             'pack_format':pack_format,
             'colnames': names,
             'dtypes': dtypes,
             'N': N,
         }, rawformat)

    # Make the storage
    fname = os.path.join(directory, "__store___")
    f = open(fname, 'wb')
    f.seek(rowbytes*N)
    f.write(b'\0')
    f.close()

    # return the store
    return RawStore(directory, Mode.WRITE)
    
    
