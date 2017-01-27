"""
This is a canned example for generating descriptors
Please modify as necessary
"""
from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors as rd
import numpy,sys
from .DescriptorGenerator import DescriptorGenerator
import logging
from rdkit.Novartis import DescriptorEngine
from rdkit import Chem

DEFAULT_DESCRIPTORS="rdkit:Chi2v,moka:fractionIonized(pH=10.0),rdkit:Kappa2,rdkit:Kappa1,rdkit:NHOHCount,rdkit:Chi2n,rdkit:HeavyAtomCount,molvol:volume,rdkit:NumHeteroatoms,rdkit:FpDensityMorgan1,rdkit:FpDensityMorgan2,rdkit:HeavyAtomMolWt,rdkit:MinAbsEStateIndex,rdkit:VSA_EState10,rdkit:BertzCT,rdkit:Kappa3,rdkit:MinAbsPartialCharge,rdkit:SMR_VSA6,rdkit:FractionCSP3,rdkit:MinEStateIndex,rdkit:SMR_VSA5,rdkit:SMR_VSA4,rdkit:SMR_VSA7,basic:NumHDonors,rdkit:SMR_VSA1,rdkit:SMR_VSA3,rdkit:MaxAbsEStateIndex,rdkit:NumAromaticHeterocycles,rdkit:SMR_VSA9,rdkit:SMR_VSA8,rdkit:Chi0,rdkit:Chi1,rdkit:NOCount,basic:NumRotatableBonds,basic:NumAmides,rdkit:Chi0v,rdkit:NumSaturatedRings,rdkit:Chi0n,rdkit:SlogP_VSA6,rdkit:SlogP_VSA7,rdkit:SlogP_VSA4,rdkit:MolLogP,rdkit:SlogP_VSA2,rdkit:SlogP_VSA3,rdkit:SlogP_VSA1,rdkit:MaxPartialCharge,rdkit:SlogP_VSA8,rdkit:SlogP_VSA9,rdkit:VSA_EState4,rdkit:qed,rdkit:VSA_EState8,rdkit:VSA_EState9,rdkit:VSA_EState6,rdkit:VSA_EState7,rdkit:Chi3v,rdkit:VSA_EState5,rdkit:VSA_EState2,rdkit:VSA_EState3,rdkit:VSA_EState1,rdkit:MaxAbsPartialCharge,rdkit:Chi3n,rdkit:MolMR,rdkit:Ipc,rdkit:TPSA,rdkit:ExactMolWt,basic:NumHAcceptors,rdkit:RingCount,rdkit:SMR_VSA2,rdkit:NumAromaticRings,rdkit:EState_VSA10,rdkit:EState_VSA11,rdkit:NumAromaticCarbocycles,rdkit:FpDensityMorgan3,rdkit:MaxEStateIndex,rdkit:MolWt,rdkit:HallKierAlpha,rdkit:Chi1v,rdkit:NumRotatableBonds,rdkit:Chi1n,rdkit:LabuteASA,rdkit:PEOE_VSA14,rdkit:PEOE_VSA11,rdkit:PEOE_VSA10,rdkit:PEOE_VSA13,rdkit:PEOE_VSA12,rdkit:NumValenceElectrons,rdkit:NumRadicalElectrons,rdkit:NumHAcceptors,rdkit:EState_VSA4,rdkit:EState_VSA5,rdkit:EState_VSA6,rdkit:EState_VSA7,basic:FlexIndex,rdkit:EState_VSA1,rdkit:EState_VSA2,rdkit:EState_VSA3,rdkit:EState_VSA8,rdkit:EState_VSA9,rdkit:SMR_VSA10,basic:mw,rdkit:MinPartialCharge,rdkit:Chi4n,rdkit:Chi4v,basic:psa,rdkit:NumHDonors,moka:fractionIonized(pH=4.0),moka:logD(pH=12.0),rdkit:NumAliphaticCarbocycles,rdkit:NumAliphaticRings,rdkit:SlogP_VSA10,rdkit:SlogP_VSA11,rdkit:SlogP_VSA12,rdkit:BalabanJ,rdkit:NumSaturatedCarbocycles,rdkit:PEOE_VSA9,rdkit:PEOE_VSA8,rdkit:SlogP_VSA5,rdkit:NumAliphaticHeterocycles,rdkit:PEOE_VSA1,rdkit:PEOE_VSA3,rdkit:PEOE_VSA2,rdkit:PEOE_VSA5,rdkit:PEOE_VSA4,rdkit:PEOE_VSA7,rdkit:PEOE_VSA6,moka:logD(pH=7.4),rdkit:NumSaturatedHeterocycles"

class DescriptorEngineDescriptors(DescriptorGenerator):
    NAME="DescriptorEngineDescriptors"
    def __init__(self, descriptors=DEFAULT_DESCRIPTORS):
        self.descriptors = descriptors.strip().split(",")
        self.engine = DescriptorEngine.GetDescriptorCalculator(self.descriptors)
        self.columns = [(name, numpy.float64) for name in self.descriptors]
        DescriptorGenerator.__init__(self)

    def GetColumns(self):
        """Returns [(name, numpy.dtype), ...] for all columns being computed"""
        return self.columns

    @staticmethod
    def get(values, key, smiles):
        r = values.get(key, None)
        if r == None:
            print("Failed to compute %s for smiles %s"%(key, smiles), file=sys.stderr)
        return float('nan')
    
    def processMol(self, m, smiles):
        res = self.engine.calculate([m])
        values = res[0]
        result = [ self.get(values, name, smiles) for name in self.descriptors ]
        return result

    def processMols(self, mols, smiles):
        res = self.engine.calculate([mols])
        results = []
        for values in res:
            values = res[0]
            result.append([ self.get(values, name, smiles) for name in self.descriptors ])
        return results
        
    def process(self, smiles):
        """smiles
        generate descriptors from a smiles string using the specified
        properties.  

        Default is to return morgan3 folded counts clipped to 255 and
        use rdkit 2D properties.
        """
        try:
            m = Chem.MolFromSmiles(smiles)
        except:
            return None

        if m == None:
            return None

        return self.processMol(m, smiles)


DescriptorEngineDescriptors()
