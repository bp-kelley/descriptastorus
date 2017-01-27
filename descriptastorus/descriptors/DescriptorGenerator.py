from rdkit import Chem
import logging

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

    def processMol(self, m, smiles):
        """rdmol, smiles
        generate descriptors from a smiles string using the specified
        properties.  

        Default is to return morgan3 folded counts clipped to 255 and
        use rdkit 2D properties.
        """
        raise NotImplementedError

    def process(self, smiles):
        """smiles
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

        return self.processMol(mol, smiles)

    def processMols(self, mols, smiles):
        return [self.processMol(m, smile) for m, smile in zip(mols, smiles)]

    def processSmiles(self, smiles):
        mols = []
        indices = []
        goodsmiles = []
        for i,smile in enumerate(smiles):
            m = self.molFromSmiles(smile)
            if m:
                mols.append(m)
                indices.append(i)
                goodsmiles.append(smile)
                
        results = processMols(mols, goodsmiles)
        if len(indices) == len(smiles):
            return result

        else:
            all_results = [None] * len(smiles)
            for idx,result in zip(indicies, results):
                all_results[idx] = result
            return all_results
        
class Container(DescriptorGenerator):
    def __init__(self, generators):
        self.generators = generators
        columns = self.columns = []
        for g in generators:
            columns.extend(g.GetColumns())

    def GetColumns(self):
        return self.columns
    
    def processMol(self, m, smiles):
        results = []
        for g in self.generators:
            results.extend(g.processMol(m,smiles))
        return results

    def processMols(self, mols, smiles):
        results = []
        for g in self.generators:
            results.extend(g.processMols(mols,smiles))
        return results
    

def MakeGenerator( generator_names ):
    generators = []
    for name in generator_names:
        try:
            d = DescriptorGenerator.REGISTRY[name.lower()]
            generators.append(d)
        except:
            logging.exception("No DescriptorGenerator found named %s", name)
            raise
    return Container(generators)
