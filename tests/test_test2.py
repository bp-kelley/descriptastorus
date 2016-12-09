from descriptastorus import MolFileIndex
index = MolFileIndex.MakeSmilesIndex("data/test2.smi", "test2", hasHeader=False,
                                     smilesColumn=1, nameColumn=0)
assert index.N == 13
assert index.getMol(12) == 'c1ccccc1CCCCCCCCCCCC'
assert index.getName(12) == '13'

