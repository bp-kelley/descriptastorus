from __future__ import print_function
from . import DescriptaStore, MolFileIndex
from .descriptors import MakeGenerator
from rdkit.Chem import AllChem

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
                 descriptors, index_inchikey, batchsize=1000):
        self.storage = storage
        self.smilesfile = smilesfile
        self.smilesColumn = smilesColumn
        self.nameColumn = nameColumn
        self.seperator = seperator
        self.descriptors = descriptors
        self.hasHeader = hasHeader
        self.index_inchikey = index_inchikey
        self.batchsize = int(batchsize)

# ugly multiprocessing nonesense
#  this makes this really not threadsafe
props = []

def process( job ):
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
    except Exception, x:
        import traceback
        traceback.print_exc()

    return res

def processInchi( job ):
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
    except Exception, x:
        import traceback
        traceback.print_exc()

    return res

# not thread safe!
def make_store(options):
    while props:
        props.pop()
        
    props.append( MakeGenerator(options.descriptors.split(",")) )
    properties = props[0]
    # to test molecule
    
    inchiKey = options.index_inchikey
    if inchiKey and not kyotocabinet:
        print("Indexing inchikeys requires kyotocabinet, please install kyotocabinet",
              file=sys.stderr)
        return False
    
    # make the storage directory
    if os.path.exists(options.storage):
        raise IOError("Directory for descriptastorus already exists: %s"%options.storage)
    
    os.mkdir(options.storage)
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
            cabinet = kyotocabinet.DB()
            inchi = os.path.join(options.storage, "inchikey.kch")
            cabinet.open(inchi, kyotocabinet.DB.OWRITER | kyotocabinet.DB.OCREATE)

        if options.nameColumn != -1:
            name_cabinet = kyotocabinet.DB()
            name = os.path.join(options.storage, "name.kch")
            name_cabinet.open(name, kyotocabinet.DB.OWRITER | kyotocabinet.DB.OCREATE)


        num_cpus = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(num_cpus)
        print ("Number of molecules to process", numstructs)

        done = False
        count = 0
        batchsize = options.batchsize
        badColumnWarning = False
        inchies = {}
        while 1:
            lastcount = count
            joblist = []
            for cpuidx in range(num_cpus):
                jobs = []

                for i in range(count, min(count+batchsize, sm.N)):
                    jobs.append((i,sm.getMol(i)))

                if i+1 > count:
                    count = i+1
                if jobs:
                    joblist.append(jobs)
            if not joblist:
                break

            if options.index_inchikey:
                results = pool.map(processInchi, joblist)
            else:
                results = pool.map(process, joblist)

            for result in results:
                if not badColumnWarning and len(result) == 0:
                    badColumnWarning = True
                    print("WARNING: no molecules processed in batch, check the smilesColumn",
                          file=sys.stderr)
                    print("WARNING: First 10 smiles:\n",
                          file=sys.stderr)
                    print("\n".join(["%i: %s"%(i,sm.get(i)) for i in range(0, min(sm.N,10))]),
                          file=sys.stderr)

                if options.index_inchikey:
                    for i,v,inchi,key in result:
                        if v:
                            s.putRow(i, v)
                        if inchi in inchies:
                            inchies[key].append(i)
                        else:
                            inchies[key] = [i]

                        if options.nameColumn != -1:
                            name = sm.getName(i)
                            if name in name_cabinet:
                                print("WARNING: name %s duplicated at molecule %d and %d"%(
                                    name, name_cabinet[name], i))
                            else:
                                name_cabinet[name] = i

                elif options.nameColumn != -1:
                    for i,v in result:
                        if v:
                            s.putRow(i, v)
                        name = sm.getName(i)
                        if name in name_cabinet:
                            print("WARNING: name %s duplicated at molecule %s and %s"%(
                                name, name_cabinet[name], i))
                        else:
                            name_cabinet[name] = i

            print("Done with %s out of %s"%(count, sm.N), file=sys.stderr)

        if options.index_inchikey:
            print("Indexing inchies", file=sys.stderr)
            for k in sorted(inchies):
                cabinet[k] = repr(inchies[k])
    finally:
        sm.close()
        s.close()
