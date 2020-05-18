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
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANimpTIES, INCLUDING, BUT NOT
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
from rdkit import Chem
import logging, numpy, sys
import pandas as pd
import pandas_flavor as pf

# set to 0 to disable caching
MAX_CACHE = 10000
cache_hit = 0
cache_miss = 0

class DescriptorGenerator:
    REGISTRY = {}
    NAME = None
    def __init__(self):
        try:
            self.REGISTRY[self.NAME.lower()] = self
        except:
            logging.exception("DescriptorGenerator must have a NAME (self.NAME)")
            raise
        # the columns to be actually calculated
        #  GetColumns returns more columns here.
        self.columns = []
        self.cache = {}
        
    def molFromSmiles(self, smiles):
        """Prepare a smiles to a molecule"""
        return Chem.MolFromSmiles(smiles)

    def molFromMol(self, mol):
        """Do any internal preperation required from a user-mol"""
        return mol
    
    def GetColumns(self):
        """Returns [(name, numpy.dtype), ...] for all columns being computed"""
        if self.NAME:
            return [ (self.NAME + "_calculated", numpy.bool) ] + self.columns
        return self.columns

    def calculateMol(self, m, smiles, internalParsing):
        """Override me for the actual calculation"""
        raise NotImplementedError
    
    def processMol(self, m, smiles, internalParsing=False):
        """rdmol, smiles -> result
        generate descriptors from a smiles string using the specified
        properties.  

        Takes the molecule as-is.  Calling this directly requires the User
        to properly prepare the molecule

        The first value returned is always True to indicate that the
        descriptors have actually been set in the store
        """
        global cache_hit, cache_miss
        # let's keep us from exploding
        if len(self.cache) > MAX_CACHE:
            self.cache.clear()

        res,_ = self.cache.get(smiles, (None,None))
        if res is not None:
            cache_hit += 1
            return res

        cache_miss += 1
        if not internalParsing:
            m = self.molFromMol(m)

        # If None's are returned from calculations
        #  mark the row as failed and replace None's with default values
        #  for storage.
        res = self.calculateMol(m, smiles, internalParsing)
        if None in res:
            logging.error("None in res")
            columns = self.GetColumns()

            for idx,v in enumerate(res):
                if v is None:
                    if self.NAME:
                        logging.error("At least one result: %s(%s) failed: %s",
                                      self.NAME,
                                      columns[idx+1][0],
                                      smiles)
                        res[idx] = columns[idx+1][1]() # default value here
                    else:
                        logging.error("At least one result: %s failed: %s",
                                      columns[idx][0],
                                      smiles)
                        res[idx] = columns[idx][1]() # default value here

            logging.info("res %r", res)
            res.insert(0, False)
        else:
            res.insert(0, True)
        self.cache[smiles] = res, m
        return res
    
    def processMols(self, mols, smiles, internalParsing=False):
        """mols, smiles -> results
        Process the molecules.  Note that smiles
        may not actually be smiles strings, but molblocks as well
        this is used for error reporting

        if internalParsing is False, takes the molecules as-is.  Otherwise
        the molecule was prepared by the DescriptorGenerator by calling the appropriate
        translation function (i.e. molFromSmiles) (e.g. used for consistently
        ordering input for MoKa descriptors)  

        Calling this directly requires the User to properly prepare the molecules if necessary
        """
        if len(mols) != len(smiles):
            raise ValueError("Number of molecules does not match number of unparsed molecules")

        result = [self.processMol(m, smile, internalParsing)
                  for m, smile in zip(mols, smiles)]
        assert len(result) == len(mols)
        return result

    def process(self, smiles):
        """smiles string -> descriptors
        returns None for invalid smiles strings

        generate descriptors from a smiles string using the specified
        properties.  

        Default is to return morgan3 folded counts clipped to 255 and
        use rdkit 2D properties.
        """
        try:
            mol = self.molFromSmiles(smiles)
        except:
            return None

        if mol == None:
            return None

        return self.processMol(mol, smiles, internalParsing=True)

    def processSmiles(self, smiles, keep_mols=True):
        """smiles -> descriptors
        Process many smiles string and generate the descriptors"""
        mols = []
        allmols = []
        indices = []
        goodsmiles = []
        _results = []
        
        for i,smile in enumerate(smiles):
            res,m = self.cache.get(smile, (None, None))
            if res:
                _results.append((i, res))
                if keep_mols:
                    allmols.append(m)
            else:
                m = self.molFromSmiles(smile)
                if m:
                    mols.append(m)
                    indices.append(i)
                    goodsmiles.append(smile)
                if keep_mols:
                    allmols.append(m)

        # all cached
        if len(_results) == len(smiles):
            all_results = [r[1] for r in _results]
            return allmols, all_results
        
        # none cached
        elif len(_results) == 0:
            results = self.processMols(mols, goodsmiles, internalParsing=True)
            
            if len(indices) == len(smiles):
                for smile, res, m in zip(smiles, results, allmols):
                    self.cache[smile] = res, m
                    
                return mols, results

            # default values are None
            all_results = [None] * len(smiles)
            for idx,result,m in zip(indices, results, allmols):
                self.cache[smiles[idx]] = result,m
                all_results[idx] = result
            return allmols, all_results
        # some cached
        else:
            results = self.processMols(mols, goodsmiles, internalParsing=True)
            all_results = [None] * len(smiles)
            # grab cached
            for i,res in _results:
                all_results[i] = res

            # grab processed
            for idx,result,m in zip(indices, results, allmols):
                self.cache[smiles[idx]] = result,m
                all_results[idx] = result
            return allmols, all_results

    def processCtab(self, ctab):
        raise NotImplementedError
    
    def processSDF(self, sdf):
        raise NotImplementedError
        
class Container(DescriptorGenerator):
    def __init__(self, generators):
        self.generators = generators
        columns = self.columns = []
        for g in generators:
            columns.extend(g.GetColumns())
        self.cache = {}

    def processMol(self, m, smiles, internalParsing=False):
        results = []
        for g in self.generators:
            results.extend(g.processMol(m,smiles, internalParsing))
            
        return results

    def processMols(self, mols, smiles, internalParsing=False):
        results = []
        for m in mols:
            results.append([])
            
        for g in self.generators:
            for result, newresults in zip(results,
                                          g.processMols(mols,smiles,
                                                        internalParsing)):
                result.extend(newresults)
        return results
    

def MakeGenerator( generator_names ):
    """Make a descriptor generator by combining multiple generators

      :param generator_names: list of available generator names

      :result: DescriptorGenerator
    """
    if not len(generator_names):
        logging.warning("MakeGenerator called with no generator names")
        raise ValueError("MakeGenerator called with no generator names")
    generators = []
    for name in generator_names:
        try:
            d = DescriptorGenerator.REGISTRY[name.lower()]
            generators.append(d)
        except:
            logging.exception("No DescriptorGenerator found named %s\nCurrently registered descriptors:\n\t%s",
                              name,
                              "\n\t".join(sorted(DescriptorGenerator.REGISTRY.keys()))
            )
            raise
    if len(generators) > 1:
        return Container(generators)
    if len(generators):
        return generators[0]

@pf.register_dataframe_method
def create_descriptors(df: pd.DataFrame,
                       mols_column_name: str,
                       generator_names: list):
    """pyjanitor style function for using the descriptor generator

    Convert a column of smiles strings or RDKIT Mol objects into Descriptors.
    Returns a new dataframe without any of the original data. This is
    intentional, as Descriptors are usually high-dimensional
    features.

    This method does not mutate the original DataFrame.

    .. code-block:: python
        import pandas as pd
        import descriptastorus.descriptors
        df = pd.DataFrame(...)
        # For "counts" kind
        descriptors = descriptastorus.descriptors.create_descriptors(
            mols_column_name='smiles', generator_names=["Morgan3Count"])
    """
    generator = MakeGenerator(generator_names)
    mols = df[mols_column_name]
    if len(mols):
        if type(mols[0]) == str:
            _, results = generator.processSmiles(mols)
        else:
            results = generator.processMols(mols, [Chem.MolToSmiles(m) for m in mols])

    else:
        results = []
    fpdf = pd.DataFrame(results, columns=generator.GetColumns())
    fpdf.index = df.index
    return fpdf

