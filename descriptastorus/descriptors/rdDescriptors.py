#  Copyright (c) 2018-2021, Novartis Institutes for BioMedical Research Inc.
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
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
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
from . import rdkit_fixes
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors as rd
import numpy
from rdkit.DataStructs import IntSparseIntVect
from rdkit.DataStructs import ConvertToNumpyArray
from .DescriptorGenerator import DescriptorGenerator
import multiprocessing as mp
import logging

import sys
logger = logging.getLogger("descriptastorus")

def to_np(vect, nbits):
    arr = numpy.zeros((nbits, ), 'i')
    ConvertToNumpyArray(vect, arr)
    return arr

def clip_sparse(vect, nbits):
    l = [0]*nbits
    for i,v in vect.GetNonzeroElements().items():
        l[i] = min(v, 255)
    return l

class Morgan(DescriptorGenerator):
    """Computes Morgan3 bitvector counts"""
    NAME = "Morgan%s"
    def __init__(self, radius=3, nbits=2048):
        if radius == 3 and nbits == 2048:
            self.NAME = self.NAME % "3"
        else:
            self.NAME = (self.NAME%radius)+"-%s"%nbits
            
        DescriptorGenerator.__init__(self)
        # specify names and numpy types for all columns
        self.radius = radius
        self.nbits = nbits
        morgan = [("m3-%d"%d, numpy.uint8) for d in range(nbits)]
        self.columns += morgan

    def calculateMol(self, m, smiles, internalParsing=False):
        counts = to_np(rd.GetMorganFingerprintAsBitVect(m,
                                                        radius=self.radius, nBits=self.nbits), self.nbits)
        return counts        

Morgan()

class MorganCounts(DescriptorGenerator):
    """Computes Morgan3 bitvector counts"""
    NAME = "Morgan%sCounts"
    def __init__(self, radius=3, nbits=2048):
        if radius == 3 and nbits == 2048:
            self.NAME = self.NAME % "3"
        else:
            self.NAME = (self.NAME%radius)+"-%s"%nbits
            
        DescriptorGenerator.__init__(self)
        # specify names and numpy types for all columns
        self.radius = radius
        self.nbits = nbits
        morgan = [("m3-%d"%d, numpy.uint8) for d in range(nbits)]
        self.columns += morgan

    def calculateMol(self, m, smiles, internalParsing=False):
        v = rd.GetHashedMorganFingerprint(m,
                                          radius=self.radius, nBits=self.nbits)
        return clip_sparse(v, self.nbits)


MorganCounts()

class ChiralMorgan(DescriptorGenerator):
    """Computes Morgan3 bitvector counts"""
    NAME = "Morgan%sCounts"
    def __init__(self, radius=3, nbits=2048):
        if radius == 3 and nbits == 2048:
            self.NAME = self.NAME % "Chiral3"
        else:
            self.NAME = (self.NAME%("Chiral%s"%radius))+"-%s"%nbits
            
        DescriptorGenerator.__init__(self)
        # specify names and numpy types for all columns
        self.radius = radius
        self.nbits = nbits
        morgan = [("cm3-%d"%d, numpy.uint8) for d in range(nbits)]
        self.columns += morgan

    def calculateMol(self, m, smiles, internalParsing=False):
        return to_np(rd.GetMorganFingerprint(
            m, radius=self.radius, nBits=self.nbits, useChirality=True), self.nbits)

ChiralMorgan()

class ChiralMorganCounts(DescriptorGenerator):
    """Computes Morgan3 bitvector counts"""
    NAME = "Morgan%sCounts"
    def __init__(self, radius=3, nbits=2048):
        if radius == 3 and nbits == 2048:
            self.NAME = self.NAME % "Chiral3"
        else:
            self.NAME = (self.NAME%("Chiral%s"%radius))+"-%s"%nbits
            
        DescriptorGenerator.__init__(self)
        # specify names and numpy types for all columns
        self.radius = radius
        self.nbits = nbits
        morgan = [("cm3-%d"%d, numpy.uint8) for d in range(nbits)]
        self.columns += morgan

    def calculateMol(self, m, smiles, internalParsing=False):
        return clip_sparse(rd.GetHashedMorganFingerprint(
            m, radius=self.radius, nBits=self.nbits, useChirality=True),
                    self.nbits)

ChiralMorganCounts()

class FeatureMorgan(DescriptorGenerator):
    """Computes Morgan3 bitvector counts"""
    NAME = "Morgan%s"
    def __init__(self, radius=3, nbits=2048):
        if radius == 3 and nbits == 2048:
            self.NAME = self.NAME % "Feature3"
        else:
            self.NAME = (self.NAME%("Feature%s"%radius))+"-%s"%nbits
            
        DescriptorGenerator.__init__(self)
        # specify names and numpy types for all columns
        self.radius = radius
        self.nbits = nbits
        morgan = [("fm3-%d"%d, numpy.uint8) for d in range(nbits)]
        self.columns += morgan

    def calculateMol(self, m, smiles, internalParsing=False):
        return to_np(rd.GetMorganFingerprintAsBitVect(
            m, radius=self.radius, nBits=self.nbits, invariants=rd.GetFeatureInvariants(m)), self.nbits)


FeatureMorgan()


class FeatureMorganCounts(DescriptorGenerator):
    """Computes Morgan3 bitvector counts"""
    NAME = "Morgan%sCounts"
    def __init__(self, radius=3, nbits=2048):
        if radius == 3 and nbits == 2048:
            self.NAME = self.NAME % "Feature3"
        else:
            self.NAME = (self.NAME%("Feature%s"%radius))+"-%s"%nbits
            
        DescriptorGenerator.__init__(self)
        # specify names and numpy types for all columns
        self.radius = radius
        self.nbits = nbits
        morgan = [("fm3-%d"%d, numpy.uint8) for d in range(nbits)]
        self.columns += morgan

    def calculateMol(self, m, smiles, internalParsing=False):
        return clip_sparse(rd.GetHashedMorganFingerprint(
            m, radius=self.radius, nBits=self.nbits, invariants=rd.GetFeatureInvariants(m)),
                    self.nbits)

FeatureMorganCounts()

class AtomPair(DescriptorGenerator):
    """Computes AtomPairs bitvector counts"""
    NAME = "AtomPairCounts"
    def __init__(self, minPathLen=1, maxPathLen=30, nbits=2048):
        if minPathLen != 1 or maxPathLen != 30 or nbits != 2048:
            self.NAME = self.NAME + ("%s-%s-%s"%(minPathLen,maxPathLen,nbits))
            
        DescriptorGenerator.__init__(self)
        # specify names and numpy types for all columns
        self.minPathLen = minPathLen
        self.maxPathLen = maxPathLen
        self.nbits = nbits
        ap = [("AP-%d"%d, numpy.uint8) for d in range(nbits)]
        self.columns += ap

    def calculateMol(self, m, smiles, internalParsing=False):
        return to_np(rd.GetAtomPairFingerprint(m, minLength=self.minPathLen, 
                                               maxLength=self.maxPathLen, nBits=self.nbits), self.nbits)


AtomPair()

class AtomPairCounts(DescriptorGenerator):
    """Computes AtomPairs bitvector counts"""
    NAME = "AtomPairCounts"
    def __init__(self, minPathLen=1, maxPathLen=30, nbits=2048):
        if minPathLen != 1 or maxPathLen != 30 or nbits != 2048:
            self.NAME = self.NAME + ("%s-%s-%s"%(minPathLen,maxPathLen,nbits))
            
        DescriptorGenerator.__init__(self)
        # specify names and numpy types for all columns
        self.minPathLen = minPathLen
        self.maxPathLen = maxPathLen
        self.nbits = nbits
        ap = [("AP-%d"%d, numpy.uint8) for d in range(nbits)]
        self.columns += ap

    def calculateMol(self, m, smiles, internalParsing=False):
        return clip_sparse(rd.GetHashedAtomPairFingerprint(m, minLength=self.minPathLen, 
                                                           maxLength=self.maxPathLen, nBits=self.nbits),
                           self.nbits)

AtomPairCounts()

class RDKitFPBits(DescriptorGenerator):
    """Computes RDKitFp bitvector"""
    NAME = "RDKitFPBits"
    def __init__(self, minPathLen=1, maxPathLen=7, nbits=2048):
        if minPathLen != 1 or maxPathLen != 7 or nbits != 2048:
          self.NAME = self.NAME + ("%s-%s-%s"%(minPathLen,maxPathLen,nbits))
            
        DescriptorGenerator.__init__(self)
        # specify names and numpy types for all columns
        self.minPathLen = minPathLen
        self.maxPathLen = maxPathLen
        self.nbits = nbits
        ap = [("RDKFP-%d"%d, numpy.uint8) for d in range(nbits)]
        self.columns += ap

    def calculateMol(self, m, smiles, internalParsing=False):
        return to_np(Chem.RDKFingerprint(m, minPath=self.minPathLen, 
                                         maxPath=self.maxPathLen, fpSize=self.nbits),
                     self.nbits)

RDKitFPBits()


class RDKitFPUnbranched(DescriptorGenerator):
    """Computes RDKitFp bitvector"""
    NAME = "RDKitUnbranchedFPBits"
    def __init__(self, minPathLen=1, maxPathLen=7, nbits=2048):
        if minPathLen != 1 or maxPathLen != 7 or nbits != 2048:
          self.NAME = self.NAME + ("%s-%s-%s"%(minPathLen,maxPathLen,nbits))
            
        DescriptorGenerator.__init__(self)
        # specify names and numpy types for all columns
        self.minPathLen = minPathLen
        self.maxPathLen = maxPathLen
        self.nbits = nbits
        ap = [("RDKFP-%d"%d, numpy.uint8) for d in range(nbits)]
        self.columns += ap

    def calculateMol(self, m, smiles, internalParsing=False):
        return to_np(Chem.RDKFingerprint(m, minPath=self.minPathLen, branchedPaths=False,
                                         maxPath=self.maxPathLen, fpSize=self.nbits),
                    self.nbits)


RDKitFPUnbranched()


RDKIT_PROPS = {"1.0.0": ['BalabanJ', 'BertzCT', 'Chi0', 'Chi0n', 'Chi0v', 'Chi1', 'Chi1n',
                         'Chi1v', 'Chi2n', 'Chi2v', 'Chi3n', 'Chi3v', 'Chi4n', 'Chi4v',
                         'EState_VSA1', 'EState_VSA10', 'EState_VSA11', 'EState_VSA2',
                         'EState_VSA3', 'EState_VSA4', 'EState_VSA5', 'EState_VSA6',
                         'EState_VSA7', 'EState_VSA8', 'EState_VSA9', 'ExactMolWt',
                         'FpDensityMorgan1', 'FpDensityMorgan2', 'FpDensityMorgan3',
                         'FractionCSP3', 'HallKierAlpha', 'HeavyAtomCount', 'HeavyAtomMolWt',
                         'Ipc', 'Kappa1', 'Kappa2', 'Kappa3', 'LabuteASA', 'MaxAbsEStateIndex',
                         'MaxAbsPartialCharge', 'MaxEStateIndex', 'MaxPartialCharge',
                         'MinAbsEStateIndex', 'MinAbsPartialCharge', 'MinEStateIndex',
                         'MinPartialCharge', 'MolLogP', 'MolMR', 'MolWt', 'NHOHCount',
                         'NOCount', 'NumAliphaticCarbocycles', 'NumAliphaticHeterocycles',
                         'NumAliphaticRings', 'NumAromaticCarbocycles', 'NumAromaticHeterocycles',
                         'NumAromaticRings', 'NumHAcceptors', 'NumHDonors', 'NumHeteroatoms',
                         'NumRadicalElectrons', 'NumRotatableBonds', 'NumSaturatedCarbocycles',
                         'NumSaturatedHeterocycles', 'NumSaturatedRings', 'NumValenceElectrons',
                         'PEOE_VSA1', 'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13',
                         'PEOE_VSA14', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5',
                         'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8', 'PEOE_VSA9', 'RingCount',
                         'SMR_VSA1', 'SMR_VSA10', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5',
                         'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10',
                         'SlogP_VSA11', 'SlogP_VSA12', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4',
                         'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7', 'SlogP_VSA8', 'SlogP_VSA9',
                         'TPSA', 'VSA_EState1', 'VSA_EState10', 'VSA_EState2', 'VSA_EState3',
                         'VSA_EState4', 'VSA_EState5', 'VSA_EState6', 'VSA_EState7', 'VSA_EState8',
                         'VSA_EState9', 'fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert', 'fr_ArN',
                         'fr_Ar_COO', 'fr_Ar_N', 'fr_Ar_NH', 'fr_Ar_OH', 'fr_COO', 'fr_COO2',
                         'fr_C_O', 'fr_C_O_noCOO', 'fr_C_S', 'fr_HOCCN', 'fr_Imine', 'fr_NH0',
                         'fr_NH1', 'fr_NH2', 'fr_N_O', 'fr_Ndealkylation1', 'fr_Ndealkylation2',
                         'fr_Nhpyrrole', 'fr_SH', 'fr_aldehyde', 'fr_alkyl_carbamate', 'fr_alkyl_halide',
                         'fr_allylic_oxid', 'fr_amide', 'fr_amidine', 'fr_aniline', 'fr_aryl_methyl',
                         'fr_azide', 'fr_azo', 'fr_barbitur', 'fr_benzene', 'fr_benzodiazepine',
                         'fr_bicyclic', 'fr_diazo', 'fr_dihydropyridine', 'fr_epoxide', 'fr_ester',
                         'fr_ether', 'fr_furan', 'fr_guanido', 'fr_halogen', 'fr_hdrzine', 'fr_hdrzone',
                         'fr_imidazole', 'fr_imide', 'fr_isocyan', 'fr_isothiocyan', 'fr_ketone',
                         'fr_ketone_Topliss', 'fr_lactam', 'fr_lactone', 'fr_methoxy', 'fr_morpholine',
                         'fr_nitrile', 'fr_nitro', 'fr_nitro_arom', 'fr_nitro_arom_nonortho',
                         'fr_nitroso', 'fr_oxazole', 'fr_oxime', 'fr_para_hydroxylation', 'fr_phenol',
                         'fr_phenol_noOrthoHbond', 'fr_phos_acid', 'fr_phos_ester', 'fr_piperdine',
                         'fr_piperzine', 'fr_priamide', 'fr_prisulfonamd', 'fr_pyridine', 'fr_quatN',
                         'fr_sulfide', 'fr_sulfonamd', 'fr_sulfone', 'fr_term_acetylene', 'fr_tetrazole',
                         'fr_thiazole', 'fr_thiocyan', 'fr_thiophene', 'fr_unbrch_alkane', 'fr_urea', 'qed']
               }

CURRENT_VERSION = "1.0.0"

FUNCS = {name:func for name, func in Descriptors.descList}
def applyFunc(name, m):
    try:
        return FUNCS[name](m)
    except:
        logger.exception("function application failed (%s->%s)",
            name, Chem.MolToSmiles(m))
                       
        return None

class RDKit2D(DescriptorGenerator):
    """Computes all RDKit Descriptors"""
    NAME = "RDKit2D"
    def __init__(self, properties=RDKIT_PROPS[CURRENT_VERSION]):
        DescriptorGenerator.__init__(self)
        # specify names and numpy types for all columns
        if not properties:
            self.columns = [ (name, numpy.float64) for name,func in sorted(Descriptors.descList) ]
        else:
            columns = self.columns
            failed = []
            
            for name in properties:
                if name in sorted(FUNCS):
                    columns.append((name, numpy.float64))
                else:
                    logger.error("Unable to find specified property %s"%name)
                    failed.append(name)
            if failed:
                raise ValueError("%s: Failed to initialize: unable to find specified properties:\n\t%s"%(
                    self.__class__.__name__,
                    "\n\t".join(failed)))
        
    def calculateMol(self, m, smiles, internalParsing=False):
        res = [ applyFunc(name, m) for name, _ in self.columns ]
        return res
    

RDKit2D()

if mp.get_start_method() == 'spawn':
    from rdkit.Chem import AllChem
    # required for tests to run, boo!
    class RDKit2DSubset(RDKit2D):
        NAME="RDKit2DSubset"
        def __init__(self):
            RDKit2D.__init__(self, properties=[
                'ExactMolWt',
                'NumAliphaticRings', 'NumAromaticCarbocycles',
                'NumAromaticHeterocycles', 'NumAromaticRings'])
    RDKit2DSubset()

    class Canonicalize(DescriptorGenerator):
        NAME="Canonicalize"
        canonicalCount = 0
        def GetColumns(self):
            return [('count', numpy.float32)]

        def molFromSmiles(self, smiles):
            m = AllChem.MolFromSmiles(smiles)
            m.SetProp("ccount", str(len(smiles)))
            return m

        def processMol(self, m, smiles, internalParsing):
            assert internalParsing == True
            return [int(m.GetProp("ccount"))]

    Canonicalize()
    
    class RDKit2DSubsetSmall(RDKit2D):
        NAME="RDKit2DSubsetSmall"
        def __init__(self):
            RDKit2D.__init__(self, properties=[
                'ExactMolWt',])

    RDKit2DSubsetSmall()

    class NanDescriptors(DescriptorGenerator):
        NAME="NANDescriptors"
        def GetColumns(self):
            return [('a', numpy.int32),
                    ('b', numpy.float32),
                    ('c', numpy.float64),
                    ('d', numpy.uint8)]

        def processMol(self, m, smiles, internalParsing=False):
            return [None]*4

    NanDescriptors()

    # this should store a calculated flag error
    class NanDescriptorsWithCalcFlags(DescriptorGenerator):
        NAME="NANDescriptorsWithCalcFlags"
        def __init__(self):
            DescriptorGenerator.__init__(self)
            self.columns =[('a', numpy.int32),
                           ('b', numpy.float32),
                           ('c', numpy.float64),
                           ('d', numpy.uint8)]

        def calculateMol(self, m, smiles, internalParsing):
            return [None]*4

    NanDescriptorsWithCalcFlags()
    
    class NanDescriptors_MS(DescriptorGenerator):
        NAME="NANDescriptors_MS"
        def GetColumns(self):
            return [('a', numpy.float32),
                    ('b', numpy.float32),
                    ('c', numpy.float32),
                    ('d', numpy.float32)]

        def processMol(self, m, smiles, internalParsing=False):
            return [float('nan')]*4
    NanDescriptors_MS()
