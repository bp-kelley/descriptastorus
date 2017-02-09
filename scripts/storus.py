#!/usr/bin/env python
"""storus
see Description below
"""
from __future__ import print_function
from descriptastorus import make_store
import argparse, logging, os, shutil, time

import sys
from rdkit import rdBase
rdBase.DisableLog("rdApp.*")

parser = argparse.ArgumentParser()
parser.add_argument("smilesfile",
                    help="file containing smiles strings")
parser.add_argument("storage",
                    help="directory in which to store the descriptors")

parser.add_argument("--descriptors", default="Morgan3Counts,RDKit2D")
parser.add_argument("--hasHeader", action="store_true",
                    help="Indicate whether the smiles file has a header row")

parser.add_argument("--index-inchikey", action="store_true",
                    help="Optionally index the descriptors with inchi keys")

#parser.add_argument("--index-smiles", action="store_true",
#                    help="Indicate whether the smiles file has a header column")

parser.add_argument("--smilesColumn", default=0,
                    help="Row index (or header name if the file has a header) for the smiles column")
parser.add_argument("--nameColumn", default=None,
                    help="Row index (or header name if the file has a header) for the name column")

parser.add_argument("--seperator", default="\t",
                    help="Row index (or header name if the file has a header) for the name column")

parser.add_argument("--batchsize", default=1000, type=int,
                    help="Batchsize for chunking up the data for processing")

parser.add_argument("--numprocs", default=-1, type=int,
                    help="Number of processers to use (-1 means use all available")

parser.add_argument("--verbose",  action="store_true",
                    help="Verbose logging")

opts = parser.parse_args()
if opts.verbose:
    logging.getLogger().setLevel(logging.INFO)


make_store.make_store(make_store.MakeStorageOptions(**vars(opts)))


        
            
            
