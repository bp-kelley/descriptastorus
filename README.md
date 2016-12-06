DescriptaStorus
===============

The descriptastorus provides fast random access to rows of properties suitable for
machine learning.

An example data set is located at:

 /db/cix/descriptastorus/magma
 
It contains (for all of magma) the folded morgan counts at radius=3 and nBits=2048
and all of the RDKit 2D descriptors implemented at the C++ layer as well as
the inchiKey and NVP indices.

[n.b.] kyotocabinet is required to read the inchiKey and name indices
  This should be installed in your environment
   (conda install kyotocabinet should do the trick)



Usage
=====

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
for key in inchiKey:
    rows.extend( d.lookupInchiKey(key) )

rows.sort()
for row in descriptors:
    descriptors = d.descriptors().get(row)
    ...
```

    
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