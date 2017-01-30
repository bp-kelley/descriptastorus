from __future__ import print_function
from rdkit import Chem
import logging, sys

class DescriptorGenerator:
    REGISTRY = {}
    def __init__(self):
        try:
            self.REGISTRY[self.NAME.lower()] = self
        except:
            logging.exception("DescriptorGenerator must have a NAME (self.NAME)")
            raise

    def molFromSmiles(self, smiles):
        return Chem.MolFromSmiles(smiles)

    def GetColumns(self):
        """Returns [(name, numpy.dtype), ...] for all columns being computed"""
        raise NotImplementedError

    def processMol(self, m, smiles, internalParsing=False):
        """rdmol, smiles -> result
        generate descriptors from a smiles string using the specified
        properties.  

        Takes the molecule as-is.  Calling this directly requires the User
        to properly prepare the molecule
        """
        raise NotImplementedError
    
    def processMols(self, mols, smiles, internalParsing=False):
        """mols, smiles -> results
        Process the molecules.  Note that smiles
        may not actually be smiles strings, but molblocks as well
        this is used for error reporting

        if internalParsing is False, takes the molecules as-is.  Otherwise
        the molecule was prepared by the DescriptorGenerator by calling the appropriate
        translation function (i.e. molFromSmiels) (i.e. consistently
        ordering input for MoKa descriptors)  

        Calling this directly requires the User
        to properly prepare the molecules if necessary
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

    def processSmiles(self, smiles):
        """smiles -> descriptors
        Process many smiles string and generate the descriptors"""
        mols = []
        indices = []
        goodsmiles = []
        for i,smile in enumerate(smiles):
            m = self.molFromSmiles(smile)
            if m:
                mols.append(m)
                indices.append(i)
                goodsmiles.append(smile)

        results = self.processMols(mols, goodsmiles, internalParsing=True)

        if len(indices) == len(smiles):
            return mols, results

        else:
            # default values are None
            all_results = [None] * len(smiles)
            for idx,result in zip(indices, results):
                all_results[idx] = result
            return mols, all_results

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

    def GetColumns(self):
        return self.columns
    
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
    if not len(generator_names):
        logging.warning("MakeGenerator called with no generator names")
        raise ValueError("MakeGenerator called with no generator names")
    generators = []
    for name in generator_names:
        try:
            d = DescriptorGenerator.REGISTRY[name.lower()]
            generators.append(d)
        except:
            logging.exception("No DescriptorGenerator found named %s", name)
            raise
    if len(generators) > 1:
        return Container(generators)
    if len(generators):
        return generators[0]

