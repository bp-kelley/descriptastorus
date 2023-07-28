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
from . import DescriptaStore, MolFileIndex
from .descriptors import MakeGenerator
from .mode import Mode
from rdkit.Chem import AllChem
import pickle
import multiprocessing
import time, os, sys, numpy, shutil
import logging
from descriptastorus import MolFileIndex, raw
from .keyvalue import KeyValueAPI
from .DescriptaStore import get_options
import tqdm

logger = logging.getLogger("descriptastorus")

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

DEFAULT_KEYSTORE = "kyotostore"
class MakeStorageOptions:
    def __init__(self, storage, smilesfile, 
                 hasHeader, smilesColumn, nameColumn, seperator,
                 descriptors, index_inchikey, batchsize=1000, numprocs=-1, verbose=False,
                 keystore=None,
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
        self.keystore = keystore or DEFAULT_KEYSTORE
        if (kw):
            logger.warning("%s: ignoring extra keywords: %r", self.__class__.__name__, kw)

# ugly multiprocessing nonesense
#  this initialized the property lists, must be called
#  prior to making the store
PROPS = []
def init_props(descriptors):
    global PROPS
    PROPS.clear()
    PROPS.append(MakeGenerator(descriptors.split(",")) )
    return PROPS

def init_props_from_store(dbdir):
    global PROPS
    options = get_options(dbdir)
    generator = MakeGenerator(options['descriptors'].split(","))
    PROPS.clear()
    PROPS.append( generator )
    return PROPS

def is_empty(result):
    if result is None:
        return True
    elif type(result) == list:
        return len(result) == 0
    
    return result.shape == 0

def process( job ):
    if job:
        logger.debug("Running on %s jobs from index %s to %s",
                        len(job), job[0][0], job[-1][0])
    else:
        logger.warning("Empty joblist")

    res = []
    try:
        smiles = [s for _,s in job]
        _, results = PROPS[0].processSmiles(smiles)

        if len(smiles) != len(results):
            logger.error("Failed batch from index %s to %s"%(
                job[0][0], job[-1][0]))
            return []
                          
        return tuple(((index, result)
                      for (index,smiles), result in zip(job, results) if not is_empty(result)))
    
    except Exception as x:
        import traceback
        traceback.print_exc()

    return res

def processInchi( job ):
    if job:
        logger.debug("Running on %s jobs from index %s to %s",
                        len(job), job[0][0], job[-1][0])
    else:
        logger.warning("Empty joblist")

    res = []
    try:
        smiles = [s for _,s in job]
        mols, results = PROPS[0].processSmiles(smiles)
        if len(smiles) != len(results):
            logger.error("Failed batch from index %s to %s"%(
                job[0][0], job[-1][0]))
            return []
        
        for i, ((index, smiles), result) in enumerate(zip(job, results)):
            m = mols[i]
            if not is_empty(result):
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
                logger.warning("Duplicated name %s at file index %s and %s",
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
    properties = init_props(options.descriptors)[0]
    # to test molecule
    
    inchiKey = options.index_inchikey
    key_value_store = None
    if inchiKey and options.keystore:
        key_value_store = KeyValueAPI.get_store(options.keystore)
        if not key_value_store:
            logger.error("Indexing inchikeys requires %s, please install or use the dbm store", options.keystore)
            return False
    
    # make the storage directory
    if os.path.exists(options.storage):
        raise IOError("Directory for descriptastorus already exists: %s"%options.storage)

    # prepare the Pool
    if options.numprocs == -1:
        num_cpus = multiprocessing.cpu_count()
    else:
        # never use more than the maximum number
        num_cpus = min(int(options.numprocs), multiprocessing.cpu_count())
            
    pool = multiprocessing.Pool(num_cpus,
                                initializer=init_props,
                                initargs=(options.descriptors,))

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
    logger.info("Creating descriptors for %s molecules...", sm.N)

                                      
    numstructs = sm.N
    s = raw.MakeStore(properties.GetColumns(), sm.N, options.storage,
                      checkDirectoryExists=False)
    try:
        if options.index_inchikey and key_value_store:
            inchi = os.path.join(options.storage, "inchikey")
            cabinet = key_value_store()
            cabinet.open(inchi, Mode.WRITE)
        else:
            logger.warning("Not logging inchi (see --index-inchkey)")
            cabinet = None

        if options.nameColumn is not None and key_value_store:
            logger.info("Creating name store")
            name = os.path.join(options.storage, "name")
            name_cabinet = key_value_store()
            name_cabinet.open(name, Mode.WRITE)
            assert os.path.exists(name_cabinet.get_actual_filename(name)), name_cabinet.get_actual_filename(name)
        else:
            logger.warning("Not storing name lookup (see --nameColumn)")
            name_cabinet = None

        logger.info("Number of molecules to process: %s", numstructs)
        
        done = False
        count = 0
        numOutput = 0
        batchsize = options.batchsize
        badColumnWarning = False
        
        inchies = {}
        names = {}
        if 1:
        #with tqdm.tqdm(total=sm.N) as pbar:
        #    pbar.set_description("Processing")
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

                #pbar.update(count-lastcount)
                procTime = time.time() - t1

                for result in results:
                    numOutput += len(result)
                    if numOutput == 0 and not badColumnWarning and len(result) == 0:
                        badColumnWarning = True
                        logger.warning("no molecules processed in batch, check the smilesColumn")
                        logger.warning("First 10 smiles:\n")
                        logger.warning("\n".join(["%i: %s"%(i,sm.get(i)) for i in range(0, min(sm.N,10))]))


                flattened = [val for sublist in results for val in sublist]
                flattened.sort()

                t1 = time.time()
                delta = 0.0
                # flatten the results so that we store them in index order
                for result in flattened:
                    if options.index_inchikey:
                        i,v,inchi,key = result
                        if not is_empty(v):
                            try:
                                s.putRow(i, v)
                            except ValueError:
                                logger.exception("Columns: %s\nData: %r",
                                                  properties.GetColumns(),
                                                  v)
                                raise
                        if inchi in inchies:
                            inchies[key].append(i)
                        else:
                            inchies[key] = [i]

                    else:
                        i,v = result
                        if not is_empty(v):
                            s.putRow(i, v)

                storeTime = time.time() - t1
                logger.debug("Done with %s out of %s.  Processing time %0.2f store time %0.2f",
                    count, sm.N, procTime, storeTime)

        if cabinet and options.index_inchikey:
            logger.info("Indexing inchies")
            t1 = time.time()
            for k in sorted(inchies):#tqdm.tqdm(sorted(inchies)):
                cabinet.set(k, inchies[k])
            logger.info("... indexed in %2.2f seconds", (time.time()-t1))
            
        if name_cabinet:
            t1 = time.time()
            logger.info("Indexing names")
            for name in sorted(names):#tqdm.tqdm(sorted(names)):
                name_cabinet.set(name, names[name])
            logger.info("... indexed in %2.2f seconds", (time.time()-t1))
    finally:
        sm.close()
        s.close()
        pool.close()
