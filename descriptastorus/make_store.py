from __future__ import print_function
from . import DescriptaStore, MolFileIndex
from .descriptors import MakeGenerator
from rdkit.Chem import AllChem
import pickle
import multiprocessing
import time, os, sys, numpy, shutil
import logging
from descriptastorus import MolFileIndex, raw
try:
    import kyotocabinet
except:
    kyotocabinet = None
    logging.warning("No kyotocabinet available")
    raise
    
# args.storage
# args.smilesfile
# args.descriptors -> descriptors to make
# args.hasHeader -> true/false for smiles file input
# args.index_inchikey -> true/false
# args.index_smiles -> index smiles strings
# args.smilesCanon -> rdkit, avalon
# args.smilesColumn
# args.nameColumn
# args.seperator

class MakeStorageOptions:
    def __init__(self, storage, smilesfile, 
                 hasHeader, smilesColumn, nameColumn, seperator,
                 descriptors, index_inchikey, batchsize=1000, numprocs=-1, verbose=False,
                 **kw):
        self.storage = storage
        self.smilesfile = smilesfile
        self.smilesColumn = smilesColumn
        self.nameColumn = nameColumn
        self.seperator = seperator
        self.descriptors = descriptors
        self.hasHeader = hasHeader
        self.index_inchikey = index_inchikey
        self.batchsize = int(batchsize)
        self.numprocs = numprocs
        self.verbose = verbose
        if (kw):
            logging.warning("%s: ignoring extra keywords: %r", self.__class__.__name__, kw)

# ugly multiprocessing nonesense
#  this makes this really not threadsafe
props = []

def process( job ):
    if job:
        logging.debug("Running on %s jobs from index %s to %s",
                        len(job), job[0][0], job[-1][0])
    else:
        logging.warning("Empty joblist")

    res = []
    try:
        smiles = [s for _,s in job]
        _, results = props[0].processSmiles(smiles)
        if len(smiles) != len(results):
            logging.error("Failed batch from index %s to %s"%(
                job[0][0], job[-1][0]))
            return []
                          

        return tuple(((index, result)
                      for (index,smiles), result in zip(job, results) if result))
    except Exception as x:
        import traceback
        traceback.print_exc()

    return res

def processInchi( job ):
    if job:
        logging.debug("Running on %s jobs from index %s to %s",
                        len(job), job[0][0], job[-1][0])
    else:
        logging.warning("Empty joblist")

    res = []
    try:
        smiles = [s for _,s in job]
        mols, results = props[0].processSmiles(smiles)
        if len(smiles) != len(results):
            logging.error("Failed batch from index %s to %s"%(
                job[0][0], job[-1][0]))
            return []
        
        for i, ((index, smiles), result) in enumerate(zip(job, results)):
            m = mols[i]
            if result:
                inchi = AllChem.MolToInchi(m)
                key = AllChem.InchiToInchiKey(inchi)
                res.append((index, result, inchi, key))

        return res
    except Exception as x:
        import traceback
        traceback.print_exc()

    return res

def getJobsAndNames(molindex, options, start, end, batchsize, nprocs, names):
    jobs = []
    for i in range(nprocs):
        jobs.append([])

    last = min(end, start+batchsize*nprocs)
    for i in range(start, last):
        moldata, name = molindex.get(i)
        if name in names:
            if options.hasHeader:
                offset = 1
            else:
                offset = 0
                logging.warning("Duplicated name %s at file index %s and %s",
                                name, names[name]+offset, i+offset)
        names[name] = i

        jobs[ i%nprocs ].append((i,moldata))
    # remove empty jobs
    jobs = [ job for job in jobs if job ]        
    return jobs, last

def getJobs(molindex, options, start, end, batchsize, nprocs):
    jobs = []
    for i in range(nprocs):
        jobs.append([])

    last = min(end, start+batchsize*nprocs)
    for i in range(start, last):
        moldata = molindex.get(i)
        jobs[ i%nprocs ].append((i,moldata))
    # remove empty jobs        
    jobs = [ job for job in jobs if job ]
    return jobs, last

# not thread safe!
def make_store(options):
    while props:
        props.pop()
        
    props.append( MakeGenerator(options.descriptors.split(",")) )
    properties = props[0]
    # to test molecule
    
    inchiKey = options.index_inchikey
    if inchiKey and not kyotocabinet:
        logging.warning("Indexing inchikeys requires kyotocabinet, please install kyotocabinet")
        return False
    
    # make the storage directory
    if os.path.exists(options.storage):
        raise IOError("Directory for descriptastorus already exists: %s"%options.storage)
    
    os.mkdir(options.storage)
    with open(os.path.join(options.storage, "__options__"), 'wb') as f:
        pickle.dump(vars(options), f)

    # index the molfile
    indexdir = os.path.join(options.storage, "__molindex__")

    sm = MolFileIndex.MakeSmilesIndex(options.smilesfile, indexdir,
                                      sep=options.seperator,
                                      hasHeader = options.hasHeader,
                                      smilesColumn = options.smilesColumn,
                                      nameColumn = options.nameColumn)
    print("Creating descriptors for %s molecules..."%sm.N)

                                      
    numstructs = sm.N
    s = raw.MakeStore(properties.GetColumns(), sm.N, options.storage,
                      checkDirectoryExists=False)
    try:
        if options.index_inchikey:
            logging.info("Creating inchi store")
            cabinet = kyotocabinet.DB()
            inchi = os.path.join(options.storage, "inchikey.kch")
            cabinet.open(inchi, kyotocabinet.DB.OWRITER | kyotocabinet.DB.OCREATE)

        if options.nameColumn is not None:
            logging.info("Creating name store")
            name_cabinet = kyotocabinet.DB()
            name = os.path.join(options.storage, "name.kch")
            name_cabinet.open(name, kyotocabinet.DB.OWRITER | kyotocabinet.DB.OCREATE)


        if options.numprocs == -1:
            num_cpus = multiprocessing.cpu_count()
        else:
            # never use more than the maximum number
            options.numprocs = min(int(options.numprocs), multiprocessing.cpu_count())
            
        pool = multiprocessing.Pool(num_cpus)
        print ("Number of molecules to process", numstructs)

        done = False
        count = 0
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
                break

            
            t1 = time.time()
            if options.index_inchikey:
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
                if options.index_inchikey:
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

        if options.index_inchikey:
            logging.info("Indexing inchies")
            t1 = time.time()
            for k in sorted(inchies):
                cabinet[k] = repr(inchies[k])
            logging.info("... indexed in %2.2f seconds", (time.time()-t1))
            
        if names:
            t1 = time.time()
            logging.info("Indexing names")
            for name in sorted(names):
                name_cabinet[name] = names[name]
            logging.info("... indexed in %2.2f seconds", (time.time()-t1))
    finally:
        sm.close()
        s.close()
