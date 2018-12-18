#!/usr/bin/env python
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

"""storus
see Description below
"""
from __future__ import print_function
from descriptastorus import DescriptaStore
import argparse, csv, logging, os, shutil, time

import sys
from rdkit import rdBase
rdBase.DisableLog("rdApp.*")

parser = argparse.ArgumentParser()

parser.add_argument("storage",
                    help="descriptastore storage")

parser.add_argument("--namefile", default="",
                    help="compounds names to export to csv, one per line")

parser.add_argument("--keep-calculatedflags", action="store_true",
                    help="With this switch enabled, the calculated flags for each desriptor subset will be kept\n(note - this turns on the --keep-missing-descriptors flag")

parser.add_argument("--keep-missing-descriptors", action="store_true",
                    help="Outputs rows with missing descriptor calculations")

parser.add_argument("--seperator", default=",",
                    help="If , outputs csv, if '\t' outputs tab seperated file")

parser.add_argument("--output-name", action="store_true", default=True,
                    help="outputs the molecule name (if available)")

# to do?
#parser.add_argument("--output-smiles", action="store_true",
#                    help="outputs the smiles column")
#parser.add_argument("--output-inchi", action="store_true",
#                    help="outputs an inchi column (if available)")

parser.add_argument("--verbose",  action="store_true",
                    help="Verbose logging")

opts = parser.parse_args()

store = DescriptaStore(opts.storage)
if opts.seperator == ",":
    writer = csv.writer(sys.stdout)
elif opts.seperator == "\t":
    writer = csv.writer(sys.stdout, dialiect=csv.excel_tab)

if opts.keep_missing_descriptors:
    opts.keep_calculatedflags = True

if opts.output_name:
    writer.writerow(['Name'] + store.getDescriptorNames(opts.keep_calculatedflags))
else:        
    writer.writerow(store.getDescriptorNames(opts.keep_calculatedflags))

if not opts.namefile:
    indices = range(len(store))
else:
    def idx_getter(namefile):
        for line in open(namefile):
            try:
                idx = store.lookupName(line.strip())
                yield idx
            except:
                logging.warning("Name %r not found in store", line.strip())
                continue
    indices = idx_getter(opts.namefile)
    
for idx in indices:
    descriptors = store.getDescriptors(idx, opts.keep_calculatedflags)
    moldata = store.index.get(idx)
    if opts.output_name:
        if len(moldata) == 2:
            _, name = moldata
        else:
            name = str(i)

        row = [name] + map(repr, descriptors)
    else:
        row = map(repr, descriptors)

    writer.writerow(row)
            
