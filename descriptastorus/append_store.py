from __future__ import print_function
from . import DescriptaStore, MolFileIndex, Mode
from .descriptors import MakeGenerator
from rdkit.Chem import AllChem
import pickle
import time, os, sys, numpy, shutil
import logging
from make_store import process, processInchi, getJobsAndNames, getJobs, props
import multiprocessing

from descriptastorus import MolFileIndex, raw
try:
    import kyotocabinet
except:
    kyotocabinet = None
    logging.warning("No kyotocabinet available")
    raise
    
class AppendStorageOptions:
    def __init__(self, storage, smilesfile, 
                 hasHeader, smilesColumn, nameColumn, seperator,
                 descriptors, index_inchikey, batchsize=1000, numprocs=-1, verbose=False,
                 **kw):
        self.storage = storage
        self.smilesfile = smilesfile
        self.smilesColumn = smilesColumn
        self.nameColumn = nameColumn
        self.seperator = seperator
        self.hasHeader = hasHeader
        self.batchsize = int(batchsize)
        self.numprocs = numprocs
        self.verbose = verbose

def check(store, opts ,member_name):
    """store, opts, member_name -> None or ValueError
    checks that the value set to member_name is the same int the storage options and
    the requested options"""
    store_v = store[member_name]
    opts_v = getattr(opts, member_name)
    if store_v != opts_v:
        raise ValueError("Requsted option %r does not match original storage. (%r != %r)"%(
            member_name, store_v, opts_v))
    

def append_smiles_file(src, dest, hasHeader):
    # check to see if we end with a new line
    endsInNewLine = False
    sz = os.path.getsize(dest)

    with open(dest) as f:
        f.seek(sz - 1)
        endsInNewLine = f.read(1) == "\n"

    # append our new data to the original
    # in the cheapest way possible...
    with open(dest, 'a') as f:
        if not endsInNewLine:
            f.write("\n")
        for i,line in enumerate(open(src)):
            if i == 0:
                if not hasHeader:
                    f.write(line)
            else:
                f.write(line)
    
# not thread safe!
def append_store(options):
    # ugly
    while props:
        props.pop()
        
    # to test molecule
    
    # make the storage directory
    if not os.path.exists(options.storage):
        raise IOError("Directory for descriptastorus does not exist: %s"%options.storage)
    

    with open(os.path.join(options.storage, "__options__")) as f:
        storageOptions = pickle.load(f)
    dbdir = options.storage
    d = DescriptaStore(dbdir, mode=Mode.READONLY)
    props.append( d.getDescriptorCalculator() )
    properties = props[0]

    if d.inchikey and not kyotocabinet:
        logging.warning("Indexing inchikeys requires kyotocabinet, please install kyotocabinet")
        return False
    
    d.close()
    # checks basic compatibility (only smiles column and name and seperator)
    check(storageOptions, options, "smilesColumn")
    check(storageOptions, options, "nameColumn")
    check(storageOptions, options, "seperator")

    # reindex with the new data
    indexdir = os.path.join(options.storage, "__molindex__")
    molindex = MolFileIndex.MolFileIndex(indexdir)
    
    origN = molindex.N
    
    orig_filename = molindex.filename
    molindex.close()
    molindex = None

    append_smiles_file(src=options.smilesfile,
                       dest=orig_filename,
                       hasHeader=options.hasHeader)
                
        
    sm = MolFileIndex.MakeSmilesIndex(orig_filename, indexdir,
                                      sep=options.seperator,
                                      hasHeader = options.hasHeader,
                                      smilesColumn = options.smilesColumn,
                                      nameColumn = options.nameColumn,
                                      reIndex=True)
    assert sm.N > origN
    smN = sm.N
    sm.close()
    sm = None

    d = DescriptaStore(dbdir, mode=Mode.APPEND)
    assert d.index.N == smN
    print("Appending descriptors for molecules %s to %s..."%(
        origN,
        d.index.N))
    d.db.appendBlankRows( d.index.N - origN)
    assert d.db.N == d.index.N
    
    start = origN
    end = smN
    numstructs = smN

    try:
        sm = d.index
        s = d.db
        if d.inchikey:
          cabinet = d.inchikey

        if d.name:
            name_cabinet = d.name

        if options.numprocs == -1:
            num_cpus = multiprocessing.cpu_count()
        else:
            # never use more than the maximum number
            options.numprocs = min(int(options.numprocs), multiprocessing.cpu_count())
            
        pool = multiprocessing.Pool(num_cpus)
        print ("Number of molecules to process", numstructs)

        done = False
        count = start
        batchsize = options.batchsize
        badColumnWarning = False
        inchies = {}
        names = {}
        while 1:
            lastcount = count
            if options.nameColumn is not None:
                joblist, count = getJobsAndNames(sm, options, count, numstructs, batchsize, num_cpus, names)
            else:
                joblist, count = getJobs(sm, options, count, numstructs, batchsize, num_cpus)
                    
            if not joblist:
                logging.debug("Stopping")
                break

            t1 = time.time()
            if d.inchikey:
                results = pool.map(processInchi, joblist)
            else:
                results = pool.map(process, joblist)

                
            procTime = time.time() - t1
            
            for result in results:
                if not badColumnWarning and len(result) == 0:
                    badColumnWarning = True
                    logging.warning("no molecules processed in batch, check the smilesColumn")
                    logging.warning("First 10 smiles:\n")
                    logging.warning("\n".join(["%i: %s"%(i,sm.get(i)) for i in range(0, min(sm.N,10))]))

                
            flattened = [val for sublist in results for val in sublist]
            flattened.sort()

            t1 = time.time()
            delta = 0.0
            # flatten the results so that we store them in index order
            for result in flattened:
                if d.inchikey:
                    i,v,inchi,key = result
                    if v:
                        try:
                            s.putRow(i, v)
                        except ValueError:
                            logging.exception("Columns: %s\nData: %r",
                                              properties.GetColumns(),
                                              v)
                            raise
                    if inchi in inchies:
                        inchies[key].append(i)
                    else:
                        inchies[key] = [i]
                elif options.nameColumn is not None:
                    i,v = result
                    if v:
                        s.putRow(i, v)
                            
            storeTime = time.time() - t1
            logging.info("Done with %s out of %s.  Processing time %0.2f store time %0.2f",
                count, sm.N, procTime, storeTime)

        if d.inchikey:
            t1 = time.time()
            for k in sorted(inchies):
                if k in cabinet:
                    l = eval(cabinet[k])
                    l += inchies[k]
                    cabinet[k] = repr(l)
                else:
                    cabinet[k] = repr(inchies[k])
            logging.info("... indexed in %2.2f seconds", (time.time()-t1))
            
        if names:
            t1 = time.time()
            for name in sorted(names):
                if name in name_cabinet:
                    logging.error("Name %s already exists in database,"
                                  " keeping idx %s (duplicate idx is %s)",
                                  name, name_cabinet[name], names[name])
                else:
                    name_cabinet[name] = names[name]
            logging.info("... indexed in %2.2f seconds", (time.time()-t1))
    finally:
        d.close()
