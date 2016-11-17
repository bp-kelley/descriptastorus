class DescriptorGenerator:
    def GetColumns(self):
        """Returns [(name, numpy.dtype), ...] for all columns being computed"""
        raise NotImplementedError

    def process(self, smiles):
        """smiles
        generate descriptors from a smiles string using the specified
        properties.  

        Default is to return morgan3 folded counts clipped to 255 and
        use rdkit 2D properties.
        """
        raise NotImplementedError
