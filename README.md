DescriptaStorus
===============

The descriptastorus provides (1) fast random access to rows of properties suitable for
machine learning and (2) fast random access to indexed molecule files

An example data set is located at:

 /db/cix/descriptastorus/magma
 
It contains (for all of magma) the folded morgan counts at radius=3 and nBits=2048
and all of the RDKit 2D descriptors implemented at the C++ layer as well as
the inchiKey and NVP indices.

[n.b.] kyotocabinet is required to read the inchiKey and name indices
  This should be installed in your environment.

There are three basic ways to use DescriptaStorus:

  1. Make a use MolFileIndex
  2. Make a DescriptaStore using a script
  3. Using a DescriptaStore to access properties

Installing
==========

```
git clone https://bitbucket.org/novartisnibr/rdkit-descriptastorus.git
cd rdkit-descriptastorus
python setup.py install
```

Note:  scripts are not currently installed, they should be used from
this directory.

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
['CC(C)(O)c1ccc(nc1)c4ccc3C=CN(Cc2ccc(F)cc2)c3c4', 'NVP-LEI449']
>>> idx.getName(1000)
'NVP-LEI449'
>>> idx.getMol(1000)
CC(C)(O)c1ccc(nc1)c4ccc3C=CN(Cc2ccc(F)cc2)c3c4'
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
c1ccccc1 NVP-1234
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

Using the magma descriptastore:

```
from descriptastorus import DescriptaStore
d = DescriptaStore("/db/cix/descriptastorus/magma")

# print out the column names
print(d.descriptors().colnames)

# this will take a while!
for moldata, descriptors in d:
    smiles, name = moldata
    descriptors # is a numpy array of data morgan3 counts + rdkit descriptors

# to iterate through only the descriptors:
for descriptors in d.descriptors():
    ...


# if indexed by inchikey
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
 >>> r = raw.MakeRawStore( columns, 2, "magma")
 >>> r.putRow(0, [45.223, 3])
```

Using an existing store
-----------------------

After creation, to open the read only storage:

```
 >>> r = raw.RawStore("magma")
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