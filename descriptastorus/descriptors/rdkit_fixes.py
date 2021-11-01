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
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors as _rdMolDescriptors

import logging
logger = logging.getLogger("descriptastorus")

if not hasattr(Descriptors, "qed"):
    # fix missing internal descriptors
    from .QED import qed
    def _FingerprintDensity(mol,func,*args,**kwargs):
        fp = func(mol, *args, **kwargs)
        if hasattr(fp,'GetNumOnBits'):
            val = fp.GetNumOnBits()
        else:
            val = fp.GetTotalVal()
        return float(val)/mol.GetNumHeavyAtoms()
    FpDensityMorgan1=lambda x:_FingerprintDensity(x,_rdMolDescriptors.GetMorganFingerprint,1)
    FpDensityMorgan2=lambda x:_FingerprintDensity(x,_rdMolDescriptors.GetMorganFingerprint,2)
    FpDensityMorgan3=lambda x:_FingerprintDensity(x,_rdMolDescriptors.GetMorganFingerprint,3)
    _descList = Descriptors._descList
    _descList.append( ("qed", qed) )
    _descList.append( ("FpDensityMorgan1", FpDensityMorgan1) )
    _descList.append( ("FpDensityMorgan2", FpDensityMorgan1) )
    _descList.append( ("FpDensityMorgan3", FpDensityMorgan1) )
    logger.info("Added missing QED, FpDensityMorgan1, FpDensityMorgan2, FpDensityMorgan3")
