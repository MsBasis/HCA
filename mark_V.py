import rdkit as rd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import MoleculeDescriptors
import pandas as pd
import numpy as np
import matplotlib as plt

molecules = {
    "Mocap": "CCOP(=S)(OCC)SCC",
    "Frumin": "CCOP(=S)(OCC)SCCSCC",
    "Rampart": "CCOP(=S)(OCC)SCSCC",
    "Thiacloprid": "C1=CC(=CC=C1C#N)C2=CSC(=N2)N(C)C",
    "Dichlorodiphenyltrichloroethane (DDT)": "C1=CC=C(C=C1)C(CCl)(CCl)C2=CC=CC=C2Cl",
    "Clopyralid": "C1=CC(=C(C=C1C(=O)O)Cl)N",
    "Acetamiprid": "CC1=NC=C(N1C2=CC=CC=C2)N(C)C",
    "Aminopyralid": "C1=CC(=C(C=C1N)Cl)C(=O)O",
    "Cadusafos": "CC(C)P(=S)SC(C)CSP(=S)(C(C)C)SC(C)C",
    "Clothianidin": "CC1=NC(=CN1C2=CC=CC=C2)N(C)C(=O)NC3=NC=CN=C3",
    "Dichlorodiphenyldichloroethane (DDD)": "C1=CC=C(C=C1)C(CCl)(CCl)C2=CC=CC=C2Cl",
    "Dichlorodiphenyldichloroethylene (DDE)": "C1=CC=C(C=C1)C(CCl)=C(Cl)C2=CC=CC=C2Cl",
    "Imidacloprid": "CC1=NC(=CN1C2=CC=CC=C2)N(C)C(=O)N3C=NC=N3",
    "Malathion": "CCOC(=O)CSP(=S)(OC)OC",
    "Methoxychlor": "COC1=CC=C(C=C1)C(CCl)(CCl)C2=CC=CC=C2OC",
    "Picloram": "C1=CC(=C(C=C1C(=O)N)Cl)N",
    "Terbufos": "CC(C)SP(=S)(OC)OC(C)C",
    "Thiamethoxam": "CC1=NC(=CN1C2=CC=CC=C2)N(C)C(=O)N3C=NC(=N3)N",
    "Triclopyr": "C1=CC(=C(C=C1Cl)Cl)C(=O)O",
    "Veltin": "CC1=CC(=O)NC(=O)N1C2=CC=CC=C2"
}

descriptorNames =  [desc[0] for desc in Descriptors._descList]
calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptorNames)

def Descriptors():
    pass







