"""
This is a canned example for generating descriptors
Please modify as necessary
"""
from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import Descriptors, MolFromSmiles, MolToSmiles
from rdkit.Chem import rdMolDescriptors as rd
try:
    from rdkit.Avalon import pyAvalonTools
except:
    pass

import numpy,sys
from .DescriptorGenerator import DescriptorGenerator
import logging
from rdkit.Novartis import DescriptorEngine
from rdkit import Chem

RDKIT_DESCRIPTORS = "basic:NumAmides,basic:NumHAcceptors,basic:NumHDonors,basic:NumRotatableBonds,basic:mw,rdkit:BalabanJ,rdkit:BertzCT,rdkit:Chi0,rdkit:Chi0n,rdkit:Chi0v,rdkit:Chi1,rdkit:Chi1n,rdkit:Chi1v,rdkit:Chi2n,rdkit:Chi2v,rdkit:Chi3n,rdkit:Chi3v,rdkit:Chi4n,rdkit:Chi4v,rdkit:EState_VSA1,rdkit:EState_VSA10,rdkit:EState_VSA11,rdkit:EState_VSA2,rdkit:EState_VSA3,rdkit:EState_VSA4,rdkit:EState_VSA5,rdkit:EState_VSA6,rdkit:EState_VSA7,rdkit:EState_VSA8,rdkit:EState_VSA9,rdkit:ExactMolWt,rdkit:FpDensityMorgan1,rdkit:FpDensityMorgan2,rdkit:FpDensityMorgan3,rdkit:FractionCSP3,rdkit:HallKierAlpha,rdkit:HeavyAtomCount,rdkit:HeavyAtomMolWt,rdkit:Ipc,rdkit:Kappa1,rdkit:Kappa2,rdkit:Kappa3,rdkit:LabuteASA,rdkit:MaxAbsEStateIndex,rdkit:MaxAbsPartialCharge,rdkit:MaxEStateIndex,rdkit:MaxPartialCharge,rdkit:MinAbsEStateIndex,rdkit:MinAbsPartialCharge,rdkit:MinEStateIndex,rdkit:MinPartialCharge,rdkit:MolLogP,rdkit:MolMR,rdkit:MolWt,rdkit:NHOHCount,rdkit:NOCount,rdkit:NumAliphaticCarbocycles,rdkit:NumAliphaticHeterocycles,rdkit:NumAliphaticRings,rdkit:NumAromaticCarbocycles,rdkit:NumAromaticHeterocycles,rdkit:NumAromaticRings,rdkit:NumHAcceptors,rdkit:NumHDonors,rdkit:NumHeteroatoms,rdkit:NumRadicalElectrons,rdkit:NumRotatableBonds,rdkit:NumSaturatedCarbocycles,rdkit:NumSaturatedHeterocycles,rdkit:NumSaturatedRings,rdkit:NumValenceElectrons,rdkit:PEOE_VSA1,rdkit:PEOE_VSA10,rdkit:PEOE_VSA11,rdkit:PEOE_VSA12,rdkit:PEOE_VSA13,rdkit:PEOE_VSA14,rdkit:PEOE_VSA2,rdkit:PEOE_VSA3,rdkit:PEOE_VSA4,rdkit:PEOE_VSA5,rdkit:PEOE_VSA6,rdkit:PEOE_VSA7,rdkit:PEOE_VSA8,rdkit:PEOE_VSA9,rdkit:RingCount,rdkit:SMR_VSA1,rdkit:SMR_VSA10,rdkit:SMR_VSA2,rdkit:SMR_VSA3,rdkit:SMR_VSA4,rdkit:SMR_VSA5,rdkit:SMR_VSA6,rdkit:SMR_VSA7,rdkit:SMR_VSA8,rdkit:SMR_VSA9,rdkit:SlogP_VSA1,rdkit:SlogP_VSA10,rdkit:SlogP_VSA11,rdkit:SlogP_VSA12,rdkit:SlogP_VSA2,rdkit:SlogP_VSA3,rdkit:SlogP_VSA4,rdkit:SlogP_VSA5,rdkit:SlogP_VSA6,rdkit:SlogP_VSA7,rdkit:SlogP_VSA8,rdkit:SlogP_VSA9,rdkit:TPSA,rdkit:VSA_EState1,rdkit:VSA_EState10,rdkit:VSA_EState2,rdkit:VSA_EState3,rdkit:VSA_EState4,rdkit:VSA_EState5,rdkit:VSA_EState6,rdkit:VSA_EState7,rdkit:VSA_EState8,rdkit:VSA_EState9,rdkit:qed"

class DescriptorEngineDescriptors(DescriptorGenerator):
    NAME="DescriptorEngineRDKitDescriptors"
    def __init__(self, descriptors=RDKIT_DESCRIPTORS):
        DescriptorGenerator.__init__(self)
        self.descriptors = descriptors.strip().split(",")
        self.engine = DescriptorEngine.GetDescriptorCalculator(self.descriptors)
        self.columns += [(name, numpy.float64) for name in self.descriptors]

    def molFromSmiles(self, smiles):
        return MolFromSmiles(pyAvalonTools.GetCanonSmiles(smiles, True))

    def molFromMol(self, mol):
        # we need to use a canonical ordering for the molecule
        #  for some versions of MoKa, this is Avalon ordering
        #  by default and for historical reasons
        return MolFromSmiles(pyAvalonTools.GetCanonSmiles(MolToSmiles(mol, True), True))
    
    @staticmethod
    def get(values, key, smiles):
        r = values.get(key, None)
        if r is None or r == 'None':
            #print("Failed to compute %s for smiles %s"%(key, smiles), file=sys.stderr)
            return float('nan')
        return r
    
    def calculateMol(self, m, smiles, internalParsing=False):
        res = self.engine.calculate([m])
        values = res.get(0, None)
        if values is None:
            return None
        result = [ self.get(values, name, smiles) for name in self.descriptors ]
        return result

    def processMols(self, mols, smiles, internalParsing=False):
        """Sending molecules in a batch to the descriptor calculator is way more efficient"""
        if not internalParsing:
            mols = [self.molFromMol(mol) for mol in mols]
        results = []
        res = self.engine.calculate(mols)
        for i in range(len(mols)):
            values = res.get(i, None)
            if values:
                results.append([self.get(values, name, smiles) for name in self.descriptors])
            else:
                results.append(None)
        return results
try:    
    DescriptorEngineDescriptors()
except Exception, e:
    logging.warning("Unable to load descriptor engine descriptors: %s", str(e))
    
MOKA_DESCRIPTORS = "moka:fractionIonized(pH=10.0),moka:fractionIonized(pH=4.0),"\
                   "moka:logD(pH=12.0),moka:logD(pH=7.4)"

class DescriptorEngineMokaDescriptors(DescriptorEngineDescriptors):
    NAME="DescriptorEngineMokaDescriptors"
    def __init__(self):
        DescriptorEngineDescriptors.__init__(self, MOKA_DESCRIPTORS)

try:
    DescriptorEngineMokaDescriptors()
except Exception, e:
    logging.warning("Unable to load descriptor engine descriptors: %s", str(e))
        
PSACALC_DESCRIPTORS="basic:psa"
class DescriptorEnginePSADescriptors(DescriptorEngineDescriptors):
    NAME="DescriptorEnginePSADescriptors"
    def __init__(self):
        DescriptorEngineDescriptors.__init__(self, PSACALC_DESCRIPTORS)

try:
    DescriptorEnginePSADescriptors()
except Exception, e:
    logging.warning("Unable to load descriptor engine descriptors: %s", str(e))
        
MOLVOL_DESCRIPTORS="molvol:volume"
class DescriptorEngineMolVolDescriptors(DescriptorEngineDescriptors):
    NAME="DescriptorEngineMolVolDescriptors"
    def __init__(self):
        DescriptorEngineDescriptors.__init__(self, MOLVOL_DESCRIPTORS)


try:
    DescriptorEngineMolVolDescriptors()
except Exception, e:
    logging.warning("Unable to load descriptor engine descriptors: %s", str(e))
