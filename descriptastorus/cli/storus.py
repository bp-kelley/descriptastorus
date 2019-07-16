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
from descriptastorus import append_store, make_store
import argparse, logging, os, shutil, time

import sys
from rdkit import rdBase
rdBase.DisableLog("rdApp.*")

parser = argparse.ArgumentParser()
parser.add_argument("smilesfile",
                    help="file containing smiles strings")
parser.add_argument("storage",
                    help="directory in which to store the descriptors")

parser.add_argument("--append", action="store_true",
                    help="Append new compounds to the smiles file (rejecting compounds with the same name)")
parser.add_argument("--append-store", action="store_true",
                    help="Append specified storage to the")

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

def main():
    opts = parser.parse_args()
    if opts.verbose:
        logging.getLogger().setLevel(logging.INFO)

    if opts.append and opts.append_store:
        logging.error("Use one of --append --append-store")
        
    if opts.append:
        append_store.append_smiles(append_store.AppendStorageOptions(**vars(opts)))
    elif opts.append_store:
        append_store.append_store(append_store.AppendStorageOptions(**vars(opts)))
    else:
        d = vars(opts)
        del d['append']
        make_store.make_store(make_store.MakeStorageOptions(**d))
