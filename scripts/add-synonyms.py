"""Adds synonyms to a descriptor store"""
import csv, sys, argparse, logging
from descriptastorus import DescriptaStore, raw

parser = argparse.ArgumentParser()
parser.add_argument("csvfile",
                    help="file containing smiles strings")
parser.add_argument("storage",
                    help="directory in which to store the descriptors")

parser.add_argument("--hasHeader", action="store_true",
                    help="Indicate whether the smiles file has a header row")
parser.add_argument("--preferredNameColumn", default=0,
                    help="Row index (or header name if the file has a header) for the preferred name column")
parser.add_argument("--nameColumn", default=1,
                    help="Row index (or header name if the file has a header) for the name (synonym) column")

parser.add_argument("--seperator", default="\t",
                    help="Row index (or header name if the file has a header) for the name column")


opts = parser.parse_args()
logging.getLogger().setLevel(logging.INFO)

logging.info("Loading store")
store = DescriptaStore(opts.storage, mode=raw.Mode.WRITE)

if opts.seperator == "\t":
    reader = csv.reader(open(opts.csvfile), dialect=csv.excel_tab)
elif opts.seperator == ",":    
    reader = csv.reader(open(opts.csvfile))
else:
    logging.error("Input file must be in csv format (tab or comma)")
    sys.exit(1)

header = None
if opts.hasHeader:
    header = reader.next()
    
        
try:
    preferred = int(opts.preferredNameColumn)
except:
    if not opts.hasHeader:
        logging.error(
            "When using text column names, hasHeader must be set to True (--hasHeader)")
        sys.exit(1)
    if opts.preferredNameColumn not in header:
        logging.error("preferred name column %s doesn't exist in the header",
                      opts.preferredNameColumn)
        sys.exit(1)
    preferred = header.index(opts.preferredNameColumn)
    
try:
    name = int(opts.nameColumn)
except:
    if not opts.hasHeader:
        logging.error(
            "When using text column names, hasHeader must be set to True (--hasHeader)")
        sys.exit(1)
    if opts.nameColumn not in header:
        logging.error("preferred name column %s doesn't exist in the header",
                      opts.nameColumn)
        sys.exit(1)
    name = header.index(opts.nameColumn)


logging.info("Updating synonyms")
namedb = store.name

class Processor:
    def __init__(self, reader, namedb):
        self.i = 0
        self.reader = reader
        self.namedb = namedb
        self.done = False
        self.count = 0
        self.chunksize=100000
        
    def processChunk(self):
        namedb = self.namedb
        reader = self.reader
        logging.info("Processing chunk from %d ( need %d entries )",
                     self.count, self.count+self.chunksize)
        jobs = []
        while len(jobs) < self.chunksize:
            try:
                row = reader.next()
                self.count += 1
            except StopIteration:
                self.done = True
                break

            
            try:
                p = row[preferred]
            except:
                logging.warning("Column not found %s: row: %r", preferred, row)
                continue

            try:
                n = row[name]
            except:
                logging.warning("Column not found %s: row: %r", name, row)
                continue

            try:
                idx = namedb[p]
            except IndexError:
                logging.error("%r not in indexed in datastore", p)
                continue
            
            # remove the sample id
            if "-NX-" in n:
                r = n.split("-NX-")[0]
                if r == p: continue
                n = r
            logging.info("mapping %s->%s == %s", n, p, idx)
            jobs.append((n,idx))
        jobs.sort()
        #for n, idx in jobs:
        #   namedb[n] = idx

p = Processor(reader, namedb)

while not p.done:
    namedb.transaction(p.processChunk)

logging.info("done Updating synonyms")    


    

