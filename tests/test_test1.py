from descriptastorus import MolFileIndex
index = MolFileIndex.MakeSmilesIndex("data/test1.smi", "test1", hasHeader=True,
                                     smilesColumn="smiles", nameColumn="name")
assert index.N == 13
assert index.getMol(12) == 'c1ccccc1CCCCCCCCCCCC'
assert index.getName(12) == '13'


