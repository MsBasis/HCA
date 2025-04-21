import rdkit as rd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

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

# - - - obliczanie danych - - - 
def euklidesowa(m):
    n = m.shape[0]
    result = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            dist = np.sqrt(np.sum((m[i] - m[j]) ** 2))
            result[i, j] = dist
            result[j, i] = dist
    return pd.DataFrame(result)

def maxDist(c1, c2, mo): #funkcja potrzebna do liczenia complete linkage
    return np.max([mo[i][j] for i in c1 for j in c2])

def minDist(c1, c2, mo): #funkcja potrzebna do liczenia single linkage
    return np.min([mo[i][j] for i in c1 for j in c2])

def HCA(mo, max_min):
    n = mo.shape[0]
    klastry = [[i] for i in range(n)]
    ID_Klastry = list(range(n))
    dis = mo.to_numpy().copy()
    np.fill_diagonal(dis, np.inf)

    links = []
    next_cluster_id = n

    while len(klastry) > 1:
        best = np.inf
        pair = None

        for i in range(len(klastry)):
            for j in range(i + 1, len(klastry)):
                d = max_min(klastry[i], klastry[j], dis)
                if d < best:
                    best = d
                    pair = (i, j)

        i, j = pair
        new_cluster = klastry[i] + klastry[j]
        links.append([ID_Klastry[i], ID_Klastry[j], best, len(new_cluster)])

        for k in sorted([i, j], reverse=True):
            del klastry[k]
            del ID_Klastry[k]

        klastry.append(new_cluster)
        ID_Klastry.append(next_cluster_id)
        next_cluster_id += 1

    return np.array(links)


# - - - Wizualizacja danych - - - 
def plot_dendrogram_Complete(linkage_matrix, labels=None):
    plt.figure(figsize=(10, 6))
    dendrogram(linkage_matrix, labels=labels)
    plt.title("Dendrogram – HCA (Complete Linkage)")
    plt.xlabel("Związki")
    plt.ylabel("Odległość")
    plt.show()
    
def plot_dendrogram_Single(linkage_matrix, labels=None):
    plt.figure(figsize=(10, 6))
    dendrogram(linkage_matrix, labels=labels)
    plt.title("Dendrogram – HCA (Single Linkage)")
    plt.xlabel("Związki")
    plt.ylabel("Odległość")
    plt.show()

mo = euklidesowa(X_scaled)

linkageComplete = HCA(mo,maxDist)
linkageSingle = HCA(mo,minDist)

plot_dendrogram_Complete(linkageComplete, labels=[str(i+1) for i in range(len(compound_names))])
plot_dendrogram_Single(linkageSingle,labels=[str(i+1) for i in range(len(compound_names))])
