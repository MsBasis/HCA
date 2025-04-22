# HCA
# Report: Molecular Descriptor Analysis and HCA using RDKit

## ✨ Project Goal

The aim of this analysis was to perform **Hierarchical Cluster Analysis (HCA)** on selected chemical compounds based on their molecular descriptors. The analysis employed the **Euclidean metric** and two cluster linkage methods: **complete linkage** and **single linkage**. The resulting dendrograms were compared and interpreted.

---

## 🔬 Input Data

The dataset consisted of 20 chemical compounds represented by **SMILES** notation. Each compound was labeled with an identifier from 1 to 20. Molecular descriptors were automatically generated using the **RDKit** library.

Out of the initial 217 descriptors, 68 that had constant values were removed. The remaining descriptors were used as input features for the analysis.

---

## 🔧 Calculation Process

### 1. Descriptor Calculation
Using the `MolecularDescriptorCalculator` object, all available RDKit descriptors were computed for each molecule. The data were saved in a `pandas.DataFrame` table format.

### 2. Data Cleaning and Scaling
Before analysis, the data were:
- cleaned from constant-value columns,
- stripped of missing values,
- standardized using the **Z-score** method.

### 3. Distance Matrix and HCA
- A **Euclidean distance** matrix was computed for all molecule pairs.
- HCA was performed using two linkage strategies:
  - **Complete linkage**: maximum distance between cluster members
  - **Single linkage**: minimum distance between cluster members

In both cases, dendrograms were generated using a custom-implemented algorithm.

---

## 📊 Analysis Results

### Dendrograms:
- **Complete linkage** produced clear clusters around 4–5 groups at a distance threshold of ~10.
- **Single linkage** formed 3 main clusters but with more elongated and less homogeneous groupings.

### Key Observations:
- Compounds 1 (Mocap), 2 (Frumin), and 3 (Rampart) clustered early, indicating high similarity.
- Compounds 5 (DDT), 11 (DDD), and 12 (DDE) formed consistent clusters in both dendrograms.
- Compound 9 (Cadusafos) joined other clusters very late, suggesting its chemical uniqueness.
- Pairing of 13 (Imidacloprid), 10 (Clothianidin), and 18 (Thiamethoxam) implies structural or functional similarity.

### Methodological Insights:
- **Complete linkage** produces tighter, more robust clusters, resistant to chaining.
- **Single linkage** merges faster but may create stretched, inconsistent clusters.

---

## 📁 Project Structure

```
.
├── main.py               # Main script with HCA computations
├── descriptors           # Output data from RDKit
├── Sprawozdanie_HCA.pdf  # PDF version of the report
└── README.md             # This report (markdown format)
```

---

## 📓 Author

Mateusz Gawin  
Student ID: 292409  
Academic year: 2024/2025

