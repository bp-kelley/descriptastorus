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
"""Raw storage of integer and floating point data.

This provides fast random access to indexed data.

r = RawStore(directory)
row = r.get(0)
row = r.get(100000)

Data can be updated using putRow when opened in append mode.
"""
from __future__ import print_function
from .mode import Mode
import pickle, numpy, os, mmap, struct, sys, shutil
import logging
import six

logger = logging.getLogger("descriptastorus")

# raw stores are little endian!

class RawStoreIter:
    def __init__(self, raw):
        self.raw = raw
        self.i = -1
        
    def next(self):
        return self.__name__(self)
    
    def __next__(self):
        self.i += 1
        if self.i >= self.raw.N:
            raise StopIteration()
        return self.raw.get(self.i)

def convert_string( v ):
    if type(v) == str:
        return str.encode(v)
    return v

if six.PY3:
    def tostr( v ):
        return v.decode("utf-8").split("\x00")[0]
else:
    def tostr( v ):
        return str(v).split("\x00")[0]
    
    
class RawStore:
    def __init__(self, directory, mode=Mode.READONLY):
        """Raw storage engine
        directory = existing directory to read or prepare the raw storage
        if N is None, an existing directory is created with N entries"""
        if not os.path.exists(directory):
            raise IOError("Directory %s for raw store does not exist"%
                          directory)
        self.directory = directory
        with open(os.path.join(directory, "__rawformat__"), 'rb') as rawformat:
            self.__dict__.update(pickle.load(rawformat))

        self.unpack = struct.Struct(self.pack_format).unpack
        
        fname = self.fname = os.path.join(directory, "__store___")
        colCache = self.colCacheDir = os.path.join(directory, "__colstore___")

        self._f = None
        self.mode = mode
        self._openfile()

    def _clearColCache(self):
        if os.path.exists(self.colCacheDir):
            shutil.rmtree(self.colCacheDir)
        
    def _openfile(self):
        fname = self.fname
        mode = self.mode
        if mode == Mode.APPEND:
            self._clearColCache()
            self._f = open(fname, 'rb+')
            access = None
        if mode == Mode.READONLY:
            self._f = open(fname, 'rb')
            access = mmap.ACCESS_READ # shared by default
        elif mode == Mode.READONCE:
            self.f = open(fname, 'rb')
        else:
            self._clearColCache()
            self._f = open(fname, 'r+b')
            access = mmap.ACCESS_WRITE

        if self._f is not None:
            self.f = mmap.mmap(self._f.fileno(), 0, access=access)
        else:
            self._f = self.f

    def close(self):
        self.f.close()
        self._f.close()

    def __len__(self):
        return self.N
    
    def __iter__(self):
        return RawStoreIter(self)

    def _resetSize(self):
        opts = None
        with open(os.path.join(self.directory, "__rawformat__"), 'rb') as rawformat:
            opts = pickle.load(rawformat)
        opts['N'] = self.N
        with open(os.path.join(self.directory, "__rawformat__"), 'wb') as rawformat:
            pickle.dump(opts, rawformat)
        self._openfile()
        
    def appendBlankRows(self, M):
        """Adds M blank rows to the store (must be opened in append mode)"""
        if self.mode != Mode.APPEND:
            raise IOError("Storage must be opened in append mode to add blank rows")
        self.close()
        _f = open(self.fname, 'rb+')        
        if M < 1:
            raise ValueError("The value of M must be positive, not %r"%M)
        self.close()
        
        logger.info("Seeking to %s",  self.rowbytes * (self.N + M))
        _f.seek(self.rowbytes * (self.N + M))
        _f.write(b'\0')
        self.N += M
        _f.close()
        logger.info("Filesize is %s", os.path.getsize(self.fname))
        self._resetSize()

    def append(self, raw):
        """Append one raw store to another"""

        if self.mode != Mode.APPEND:
            raise IOError("Storage must be opened in append mode to add blank rows")
            
        if raw.getColFormats() != self.getColFormats():
            raise ValueError("Column formats are not compatible for appending store")
        N = self.N
        M = len(raw)
        self.N += M
        
        self.close()
        with open(self.fname, 'r+b') as dst:
            dst.seek(N*self.rowbytes)
            with open(raw.fname, 'rb') as src:
                shutil.copyfileobj(src, dst)
        print("****** size is now", self.N)
        self._resetSize()
        print("****** size is now(2)", self.N)
        
        
    def getDict(self, idx):
        """{colname:value, ...} Return the row at idx as a dictionary"""
        return {name:v for name, v in zip(self.colnames, self.get(idx))}
        
    def get(self, idx):
        """Return the row at idx"""
        if idx >= self.N or idx < 0:
            raise IndexError("Index out of range %s (0 < %s)"%(
                             idx, self.N))
        offset = idx * self.rowbytes
        try:
            self.f.seek(offset,0)
        except ValueError:
            print("Could not seek to index %d at offset %f offset"%(
                idx, offset), file=sys.stderr)
            raise IndexError("out or range %d"%idx)
        
        _bytes = self.f.read(self.rowbytes)
        res = self.unpack(_bytes)#struct.unpack(self.pack_format, _bytes)
        if "s" not in self.pack_format:
            return res
        else:
            return tuple([ tostr(x)
                           if isinstance(x, (str, bytes)) else x for x in res ])

    def getEndian(self):
        if self.pack_format[0] in "@<>!=":
            return self.pack_format[0]
        return "@"
    
    def getColFormats(self):
        start = 0
        endian = self.getEndian()
        if self.pack_format[0] == endian:
            start = 1
        else:
            start = 0
        num = []
        types = []
        while start < len(self.pack_format):
            c = self.pack_format[start]
            if c in "0123456789":
                num.append(c)
            else:
                types.append( "".join(num) + c )
            start += 1

        return endian, types

    def getFormatAndBytesForColumn(self, i):
        endian, formats = self.getColFormats()
        bytes = struct.calcsize(formats[i])
        return formats[i], bytes

    def getOffsetToColumn(self,i):
        if i == 0: return 0
        if i == len(self.colnames):
            return self.rowbytes # sentinel
        if i > len(self.colnames):
            raise IndexError("Column index %s out of range", i)
        
        endian, types = self.getColFormats()
        pack_format = "%s%s"%(endian, "".join(types[:i+1]))

        # because of alignment issues, calcualte the whole rowsize
        #  from the beginnig then remove the last column format
        bytes = len(struct.pack(pack_format,
                                *[dtype(0)
                                  for dtype in self.dtypes[:i+1]]))
        bytes -= struct.calcsize(types[i])
        return bytes
        
            
    def getColByIdx(self, column, forceRead=False):
        """Return the data in the entire column (lazy generator)"""
        # figure out the column
        fn = os.path.join(self.colCacheDir, str(column))
        if not forceRead and os.path.exists(fn):
            # read column from cache
            with open(fn, 'rb') as f:
                endian, types = self.getColFormats()
                pack_format = "%s%s%s"%(endian, self.N,types[column])
                for v in struct.unpack(pack_format, f.read()):
                    yield v
                return
            
        # figure out how many bytes to skip per row
        #  this is a bug in the pack_format, it is broken for strings...
        if "s" in self.pack_format:
            raise TypeError("Cannot get columns for data stores with strings currently")

        offset = self.getOffsetToColumn(column)
        pack_format, nbytes = self.getFormatAndBytesForColumn(column)
        assert nbytes

        # compute the next entry
        # use absolute indexing (SEEK_SET)
        self.f.seek(0, os.SEEK_SET)
        # read the first row
        _bytes = self.f.read(self.rowbytes)
        res = struct.unpack(self.pack_format, _bytes)
        
        self.f.seek(offset, os.SEEK_SET)

        count = 0
        while count < self.N:
            bytes = self.f.read(nbytes)

            if bytes == '':
                raise StopIteration
            try:
                v =  struct.unpack(pack_format, bytes)[0]
                if count == 0:
                    # first column sanity check
                    assert v == res[column], "%r %r %r %r"%(v,
                                                            column,
                                                            res[column], res)
                count += 1
                    
                yield v
            except struct.error:
                if count != self.N:
                    logger.exception("Bad Did not read %s entries got %s", self.N, count)
                if len(bytes) != nbytes:
                    raise StopIteration
                raise

            try:
                self.f.seek(offset, os.SEEK_SET)
            except ValueError:
                if count != self.N:
                    logger.exception("Did not read %s entries got %s", self.N, count)
                raise StopIteration

            
            offset += self.rowbytes

    def cacheColumns(self):
        self._clearColCache()
        if not os.path.exists(self.colCacheDir):
            os.mkdir(self.colCacheDir)
        for idx,colname in enumerate(self.colnames):
            fn = os.path.join(self.colCacheDir, str(idx))

            endian, types = self.getColFormats()
            pack_format = "%s%s%s"%(endian, self.N, types[idx])
            column = list(self.getColByIdx(idx, forceRead=True))
            with open(fn, 'wb') as f:
                f.write( struct.pack( pack_format, *column ) )


    def writeColIdx(self, column, data):
        # figure out how many bytes to skip per row
        skip_format = self.pack_format[:column]
        offset = len(struct.pack(skip_format,
                                 *[dtype(0)
                                   for dtype in self.dtypes[:column]]))


        # figure out the bytes to read via the format for the row
        pack_format = self.pack_format[column]
        nbytes = len(struct.pack(pack_format, self.dtypes[column](0)))

        # seek to the first column

        for v in data:
            try:
                self.f.seek(offset,0)
            except ValueError:
                break
            
            bytes = struct.pack(pack_format, v)
            self.f.write(bytes)
            offset += self.rowbytes
            
        
    def getCol(self, column_name):
        """Return the column addressed by column_name
        throws IndexError if the column doesn't exist."""
        idx = self.colnames.index(column_name)
        return self.getColByIdx(idx)

    def writeColByIdx(self, column, data):
        # figure out how many bytes to skip per row
        skip_format = self.pack_format[:column]
        offset = len(struct.pack(skip_format,
                                 *[dtype(0)
                                   for dtype in self.dtypes[:column]]))


        # figure out the bytes to read via the format for the row
        pack_format = self.pack_format[column]
        nbytes = len(struct.pack(pack_format, self.dtypes[column](0)))

        # seek to the first column

        for v in data:
            try:
                self.f.seek(offset,0)
            except ValueError:
                break
            
            bytes = struct.pack(pack_format, v)
            self.f.write(bytes)
            offset += self.rowbytes

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
        try:
            v = [dtype(v) for dtype,v in zip(self.dtypes, row)]
        except TypeError as e:
            # we might have None's here
            message = [str(e)]
            for i,(dtype,v) in enumerate(zip(self.dtypes, row)):
                try: dtype(v)
                except:
                    message.append("\tFor column %r can't convert %r to %s"%(
                        self.colnames[i],
                        v,
                        dtype))
            raise TypeError("\n".join(message))
        
        offset = idx * self.rowbytes
        self.f.seek(offset,0)
        try:
            bytes = struct.pack(self.pack_format, *[convert_string(x) for x  in row])
        except struct.error:
            logger.exception("Can't write row %r\ntypes: %r\nformat: %r",
                              row,
                              self.dtypes,
                              self.pack_format)
            raise
        try:
            self.f.write(bytes)
        except Exception as e:
            logger.error("Attempting to write to offset: %s", offset)
            logger.error("Rowsize: %s", self.rowbytes)
            logger.error("Row: %s", idx)
            logger.error("Max Row: %s", self.N)
            logger.error("Filesize is %s",
                           os.path.getsize(self.fname))
            raise

    def write(self, row):
        """Writes row with no datatype checking to the end of the file.
        Do not use with putRow/getRow
        Rows must be written in the correct order
        (strings datatypes are currently an issue)
        """
        bytes = struct.pack(self.pack_format, *row)
        self.f.write(bytes)

def str_store(s):
    return str.encode(str(s))
        
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
        elif dtype == numpy.float16:
            type = "e"
            dtypes.append(float)
        elif dtype == numpy.float32:
            type = "f"
            dtypes.append(float)
        elif dtype == numpy.float64:
            type = "d"
            dtypes.append(float)
        elif dtype == bool:
            type = "?"
            dtypes.append(bool)
        elif hasattr(dtype, 'type'): # for strings
            if dtype.type == numpy.string_:
                size = dtype.itemsize
                type = "%ss"%size
                dtypes.append(str_store)
        else:
            raise ValueError("Unhandled numpy type %s"%dtype)

        types.append(type)
                     

    # default to unaligned bytes since we are saving to disk
    pack_format = "<" + "".join(types)
    #opack_format = "".join(types)    
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
    with open(fname, 'wb') as f:
        f.seek(rowbytes*N)
        f.write(b'\0')

    # return the store
    return RawStore(directory, Mode.WRITE)
    
    
    
