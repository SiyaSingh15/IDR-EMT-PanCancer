import pandas as pd
import numpy as np
import os

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

# Load expression matrix
print("Loading expression matrix...")
expr = pd.read_csv("data/expression/emt_expression_matrix.tsv", sep="\t")
print(f"Shape: {expr.shape}")
print(f"Samples: {expr['sample_id'].nunique():,}")
print(f"Cancer types: {expr['cancer_type'].nunique()}")
print(f"Columns (first 10): {list(expr.columns[:10])}")

# Set sample_id as index
expr = expr.set_index("sample_id")

# Define gene signatures

# MSigDB Hallmark EMT — mesenchymal genes (high = mesenchymal)
MESEN_GENES = [
    "ACTA2","CDH2","COL1A1","COL1A2","FN1","CTGF","FERMT2","FGF2",
    "FLNA","IGFBP3","IL6","ITGA5","ITGB1","ITGB3","LAMA3","LAMC2",
    "MMP2","MMP3","MMP14","POSTN","RHOC","SPARC","TAGLN","TGFB1",
    "TGFB3","TGFBI","THBS1","THBS2","TIMP1","TIMP3","TNC","TPM1",
    "TWIST1","TWIST2","VIM","WNT5A","ZEB1","ZEB2","SNAI1","SNAI2",
    "AXL","MET","S100A4","YAP1","VEGFA","HIF1A","STAT3","MYC"
]

# Epithelial genes (high = epithelial)
EPITH_GENES = [
    "CDH1","EPCAM","KRT18","KRT19","KRT8","ESRP1","ESRP2",
    "GRHL2","OVOL1","OVOL2","CD24","CTNNB1"
]

# Hybrid E/M marker genes (Jolly EMP signature)
# High in hybrid state — neither fully E nor fully M
HYBRID_GENES = ["CD44","EGFR","MET","AXL","CXCR4","ITGB3","ITGB1",
                "CD24","CTNNB1","MYC"]

# Compute scores
def mean_signature_score(df, genes, name):
    """
    Simple mean expression score across signature genes.
    Robust, interpretable, and standard for bulk RNA-seq.
    Genes not in df are skipped with a warning.
    """
    available = [g for g in genes if g in df.columns]
    missing = [g for g in genes if g not in df.columns]
    if missing:
        print(f"  [{name}] Missing genes: {missing}")
    score = df[available].mean(axis=1)
    return score

print("\nComputing EMT scores...")

# Mesenchymal score (high = more mesenchymal)
expr["mesen_score"] = mean_signature_score(expr, MESEN_GENES, "Mesenchymal")

# Epithelial score (high = more epithelial)
expr["epith_score"] = mean_signature_score(expr, EPITH_GENES, "Epithelial")

# EMP score = mesenchymal - epithelial
# Positive = mesenchymal, Negative = epithelial, Near zero = hybrid
expr["emp_score"] = expr["mesen_score"] - expr["epith_score"]

print("Done.")

# Classify tumours into E / Hybrid / M states
# Bottom 25% = Epithelial, Top 25% = Mesenchymal, Middle 50% = Hybrid
q25 = expr["emp_score"].quantile(0.25)
q75 = expr["emp_score"].quantile(0.75)

print(f"\nEMP score quartiles:")
print(f"  Q25 (E/Hybrid boundary):  {q25:.3f}")
print(f"  Q75 (Hybrid/M boundary):  {q75:.3f}")

def classify_emt(score):
    if score <= q25:
        return "Epithelial"
    elif score >= q75:
        return "Mesenchymal"
    else:
        return "Hybrid"

expr["emt_state"] = expr["emp_score"].apply(classify_emt)

print(f"\nEMT state distribution:")
print(expr["emt_state"].value_counts())

print(f"\nEMT state by cancer type:")
state_by_cancer = expr.groupby("cancer_type")["emt_state"].value_counts().unstack(fill_value=0)
print(state_by_cancer.to_string())

# Save EMT score matrix
score_cols = ["cancer_type", "mesen_score", "epith_score", "emp_score", "emt_state"]
scores_df = expr[score_cols].copy()
scores_df.index.name = "sample_id"

out = "results/emt_scores.tsv"
scores_df.to_csv(out, sep="\t")
print(f"\nSaved EMT scores to: {out}")
print(f"Total samples scored: {len(scores_df):,}")
