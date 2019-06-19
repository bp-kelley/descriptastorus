DescriptaStorus
===============

The descriptastorus provides 

  1. fast random access to rows of properties suitable for
machine learning and 
  2. fast random access to indexed molecule files
  3. A mechanism for generating new descriptors for new molecules
  4. A mechanism for validating that you can recreate the same storage in different software/hardware environments
  5. An easy script for making your own descriptor files from raw data.

[n.b.] kyotocabinet is required to read/write the inchiKey and name indices
  This should be installed in your environment.

There are three basic ways to use DescriptaStorus:
  
  1. Make a DescriptaStore using a script
  2. Append new data to the store
  3. Use a DescriptaStore to access properties

Installing
==========

```
1. install rdkit
2. install scikit-learn
pip install git+https://github.com/bp-kelley/descriptastorus
```

Requirements are in the setup.py file, but essentially:

 1. python2/3
 2. rdkit
 3. [optional but highly recommended] kyotocabinet

Using RDKit descriptors
=======================
Grab a descriptor generator from the registry.

Currently registered descriptors:

	* atompaircounts
	* morgan3counts
	* morganchiral3counts
	* morganfeature3counts
	* rdkit2d
	* rdkit2dnormalized
	* rdkitfpbits

Descriptors are input as a tuple or list to the generator.

```
from descriptastorus.descriptors.DescriptorGenerator import MakeGenerator
generator = MakeGenerator(("RDKit2D",))
for name, numpy_type in generator.GetColumns():
  print("name: {} data type: {}".format(name, numpy_type))
```

The resulting columns and datatypes look like:
```
name: RDKit2D_calculated data type: <class 'bool'>
name: BalabanJ data type: <class 'numpy.float64'>
name: BertzCT data type: <class 'numpy.float64'>
name: Chi0 data type: <class 'numpy.float64'>
name: Chi0n data type: <class 'numpy.float64'>
name: Chi0v data type: <class 'numpy.float64'>
name: Chi1 data type: <class 'numpy.float64'>

```

Note: RDKit2D_calculated is just a flag for the store to indicate that the
RDKit2D features were successfully calculated.

To get combine multiple generators simply add them to the list
of desired datatypes:

```
from descriptastorus.descriptors.DescriptorGenerator import MakeGenerator
generator = MakeGenerator(("RDKit2D", "Morgan3Counts"))
smiles = "c1ccccc1"
data = generator.process(smiles)
assert data[0] is True
```

The first element is True if the molecule was successfully processed, this is used
in the descriptastor to indicate that the row is valid.

If a molecule is unsuccessfully processed, None is returned

```
data = generator.process("not a smiles")
assert data is None
```

Individual descriptor sets can also be used outside of the
generator.

```
from descriptastorus.descriptors import rdNormalizedDescriptors
from rdkit import Chem
import logging

# make the normalized descriptor generator
generator = rdNormalizedDescriptors.RDKit2DNormalized()
generator.columns # list of tuples:  (descriptor_name, numpytype) ...

# features = generator.process(smiles)
# features[0] is True/False if the smiles could be processed correcty
# features[1:] are the normalized descriptors as described in generator.columns

# example for converting a smiles string into the values
def rdkit_2d_normalized_features(smiles: str):
    # n.b. the first element is true/false if the descriptors were properly computed
    results = generator.process(smiles)[
    processed, features = results[0], results[1:]
    if processed is None:
       logging.warning("Unable to process smiles %s", smiles)
    # if processed is None, the features are are default values for the type
    return features
```

Making a DescriptaStore
=======================

see scripts/storus.py for more details:

```
usage: storus.py [-h] [--hasHeader] [--index-inchikey]
                 [--smilesColumn SMILESCOLUMN] [--nameColumn NAMECOLUMN]
                 [--seperator SEPERATOR]
                 smilesfile storage

positional arguments:
  smilesfile            file containing smiles strings
  storage               directory in which to store the descriptors

optional arguments:
  -h, --help            show this help message and exit
  --hasHeader           Indicate whether the smiles file has a header row
  --index-inchikey      Optionally index the descriptors with inchi keys
  --smilesColumn SMILESCOLUMN
                        Row index (or header name if the file has a header)
                        for the smiles column
  --nameColumn NAMECOLUMN
                        Row index (or header name if the file has a header)
                        for the name column
  --seperator SEPERATOR
                        Row index (or header name if the file has a header)
                        for the name column

```

Example:

Suppose you have a smiles file like the following:

```
SMILES STRU_ID
c1ccccc1 NAME
```

This is a whitespace seperated file with a header.  To make the standard
storage and also index the inchikey:

```
python scripts/storus.py --smilesColumn=SMILES --nameColumn=STRU_ID --hasHeader --index-inchikey \
  --seperator=" " \
  smiles.txt mysmiles-store
```

Note that smiles files are very seperator dependent.  If the smiles or name column
can't be found, it is might be because the seperator is misspecified.

The default properties created are 'Morgan3Counts,RDKit2D'.

Using a DescriptaStore
======================

Using the descriptastore (the descriptastore is a directory of files):

```
from descriptastorus import DescriptaStore
d = DescriptaStore("/db/cix/descriptastorus/store")

# print out the column names
print(d.descriptors().colnames)

# this will take a while!
for moldata, descriptors in d:
    smiles, name = moldata
    descriptors # is a numpy array of data morgan3 counts + rdkit descriptors
```

Note that the descriptors may contain status flags named as "X_Calculated" where X
is one of the descriptor sets, such as RDKit2D.

These are not returned by the iterator, or through the following api points:

```
colnames = d.getDescriptorNames()
descriptors = d.getDescriptors(index)
for moldata, descriptors in d:
  ...
```

To obtain these flags, you can either set the keepCalculatedFlags option

```
colnames = d.getDescriptorNames(keepCalculatedFlags=True)
descriptors = d.getDescriptors(keepCalculatedFlags=True)
```

or use the direct descriptor interface:

```
# to iterate through only the descriptors:
for descriptors in d.descriptors():
    ...
```

# to lookup by name (requires kyotocabinet)

```
rows = []
for name in names:
    rows.extend( d.lookupName(name) )

# sorting the rows helps with disk seeking
rows.sort()
for row in rows:
    descriptors = d.getDescriptors(row)
    ...
```

# To lookup by inchikey (requires kyotocabinet)

```
rows = []
for key in inchiKeys:
    rows.extend( d.lookupInchiKey(key) )

rows.sort()
for row in rows:
    descriptors = d.descriptors().get(row)
    smiles, name = d.molIndex().get(row)
    ...
```

Doing things yourself
=====================
    
Creating a Raw store
--------------------

The storage system is quite simple.  It is made by specifying the column names and
numpy types to store and also the number of rows to initialize.

Example:

```
 >>> from descriptastorus import raw
 >>> import numpy
 >>> columns = [('exactmw', numpy.float64), ('numRotatableBonds', numpy.int32) ...]
 >>> r = raw.MakeStore( columns, 2, "store")
 >>> r.putRow(0, [45.223, 3])
```

Using an existing store
-----------------------

After creation, to open the read only storage:

```
 >>> r = raw.RawStore("store")
```

Get the number or rows:

```
 >>> r.N
 2
```

Get the column names:

```
 >>> r.colnames
 ['exactmw', 'numRotatableBonds']
```

Extract the column:

```
>>> r.get(0)
[45.223, 3]
```

Make a MolFileIndex
===================

If the smiles file has a header

```
>>> from descriptastorus import MolFileIndex
>>> index = MolFileIndex.MakeSmilesIndex("data/test1.smi", "test1", hasHeader=True,
...                                      smilesColumn="smiles", nameColumn="name")
>>> index.N
13
>>> index.getMol(12)
'c1ccccc1CCCCCCCCCCCC'
>>> index.getName(12)
13
```

If the smiles file has no header

```
>>> from descriptastorus import MolFileIndex
>>> index = MolFileIndex.MakeSmilesIndex("data/test2.smi", "test2", hasHeader=False,
...                                      smilesColumn=1, nameColumn=0)
>>> index.N
13
>>> index.getMol(12)
'c1ccccc1CCCCCCCCCCCC'
>>> index.getName(12)
13
```

Use a MolFileIndex
==================

Using a molfile index is fairly simple:

```
>>> from descriptastorus import MolFileIndex
>>> idx = MolFileIndex("/db/cix/descriptastorus/test")
>>> idx.get(1000)
['CC(C)(O)c1ccc(nc1)c4ccc3C=CN(Cc2ccc(F)cc2)c3c4', 'XX-AAAA']
>>> idx.getName(1000)
'XX-AAAA'
>>> idx.getMol(1000)
CC(C)(O)c1ccc(nc1)c4ccc3C=CN(Cc2ccc(F)cc2)c3c4'
```

Installation
============

```
  git clone https://bitbucket.org/novartisnibr/rdkit-descriptastorus.git
  cd rdkit-descriptastorus
  python setup.py install
```


TODO:

  * fast forwards iteration (fast now, but could be faster)
  * faster append-only store creation
  * Fast molecule indexing/lookup (almost done)
  * Output to bcolz pandas backend
