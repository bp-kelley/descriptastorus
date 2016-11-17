class SmilesIndex:
    def __init__(self, filename, db,
                 smilesColumn,
                 nameColumn=-1,
                 has_header=False,
                 sep=None):
        self.db = db
        self.filename = filename
        self.has_header = has_header
        self.smilesColumn = smilesColumn
        self.nameColumn = nameColumn
        
        self.sep = sep
        # mmap?
        self.f = open(self.filename, 'r')
        
        if self.has_header:
            colnames = self.colnames = self.get(0)
            print "*"*44
            print colnames
        else:
            colnames = self.colnames = ["column_%d"%x for x in range(len(self.get(0)))]
            
        try:
            self.colidx = int(self.smilesColumn)
        except ValueError:
            self.colidx = colnames.index(self.smilesColumn)
            if self.colidx == -1:
                raise IndexError("Specified smiles column %r name not in header"%
                                 self.smilesColumn)

        if self.nameColumn != -1:
            try:
                self.nameidx = int(self.nameColumn)
            except ValueError:
                self.nameidx = colnames.index(self.nameColumn)

            if self.nameidx == -1:
                raise IndexError("Specified name column name %r not in header"%
                                 self.nameColumn)

        # get the first entry
        if self.has_header:
            row = self.get(1)
        else:
            row = self.get(0)
            
        if len(row) >= self.colidx:
            raise IndexError("Smiles Column %d greater than rowsize %s\n"
                             "Perhaps the seperator is mispecified (currently %r)"% (
                                 self.colidx,
                                 len(row),
                                 self.sep))

    def get(self, idx):
        start = self.db.get(idx)[0]
        end = self.db.get(idx+1)[0]
        self.f.seek(start,0)
        buf = self.f.read(end-start)
        print buf
        print start, end, buf
        return buf.split(self.sep)
    
    def getSmiles(self, idx):
        if self.has_header:
            idx += 1
        return self.get(idx)[self.colidx]

    def getName(self, idx):
        if self.nameidx == -1:
            raise ValueError("SmilesIndex does not have a name column")
        
        if self.has_header:
            idx += 1
        return self.get(idx)[self.nameidx]
