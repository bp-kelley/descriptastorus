from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors as _rdMolDescriptors

import logging

if not hasattr(Descriptors, "qed"):
    # fix missing internal descriptors
    from .QED import qed
    def _FingerprintDensity(mol,func,*args,**kwargs):
        fp = func(*((mol,)+args),**kwargs)
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
    logging.info("Added missing QED, FpDensityMorgan1, FpDensityMorgan2, FpDensityMorgan3")
