import rdkit as rd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


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
# - - - przygotowanie danych - - - 
descriptorNames =  [desc[0] for desc in Descriptors._descList]
calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptorNames)

def Descriptors(mols):
    results = []
    for name, smi in mols.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print("Error")
            desc = [None] * len(descriptorNames)
        else:
            desc = calculator.CalcDescriptors(mol)
        results.append([name]+list(desc))
    df = pd.DataFrame(results, columns=["Compound"]+descriptorNames)
    return df

desc_df = Descriptors(molecules)
#print(desc_df.head())
#desc_df.to_csv("Uni_Shits/HCA/descriptors",index=False)

def clean_and_scale(df):
    df_cleaned = df.copy()
    if "Compound" in df_cleaned.columns:
        compound_names = df_cleaned["Compound"]
        df_cleaned = df_cleaned.drop(columns=["Compound"])
    else:
        compound_names = pd.Series([f"Mol_{i}" for i in range(len(df_cleaned))])
    df_cleaned = df_cleaned.loc[:, df_cleaned.nunique() > 1]
    df_cleaned = df_cleaned.dropna(axis=1)
    mean = df_cleaned.mean()
    std = df_cleaned.std(ddof=0)
    X_scaled = (df_cleaned - mean) / std
    scaled_df = X_scaled.copy()
    scaled_df.insert(0, "Compound", compound_names.values)

    return scaled_df, X_scaled.to_numpy(), compound_names

scaled_df, X_scaled, compound_names = clean_and_scale(desc_df)




