import requests
import pandas as pd
import os
import time

DATA_DIR = "data/expression"
os.makedirs(DATA_DIR, exist_ok=True)

BASE = "https://www.cbioportal.org/api"

# EMT genes with their Entrez IDs
EMT_GENES = {
    "ACTA2": 59, "CDH2": 1000, "COL1A1": 1277, "COL1A2": 1278,
    "FN1": 2335, "CTGF": 1490, "CXCR4": 7852, "FERMT2": 10979,
    "FGF2": 2247, "FLNA": 2316, "IGFBP3": 3486, "IL6": 3569,
    "ITGA5": 3678, "ITGB1": 3688, "ITGB3": 3690, "LAMA3": 3909,
    "LAMC2": 3918, "MMP2": 4313, "MMP3": 4314, "MMP14": 4323,
    "POSTN": 10631, "RHOC": 389, "SPARC": 6678, "TAGLN": 6876,
    "TGFB1": 7040, "TGFB3": 7043, "TGFBI": 7045, "THBS1": 7057,
    "THBS2": 7058, "TIMP1": 7076, "TIMP3": 7078, "TNC": 3371,
    "TPM1": 7168, "TWIST1": 7291, "TWIST2": 117581, "VIM": 7431,
    "WNT5A": 7474, "ZEB1": 6935, "ZEB2": 9839, "SNAI1": 6615,
    "SNAI2": 6591, "CDH1": 999, "EPCAM": 4072, "KRT18": 3875,
    "KRT19": 3880, "KRT8": 3856, "ESRP1": 54845, "ESRP2": 80004,
    "GRHL2": 79977, "OVOL1": 5017, "OVOL2": 58495, "AXL": 558,
    "MET": 4233, "EGFR": 1956, "CD44": 960, "CD24": 100133941,
    "S100A4": 6275, "CTNNB1": 1499, "MYC": 4609, "VEGFA": 7422,
    "HIF1A": 3091, "STAT3": 6774, "YAP1": 10413,
}

STUDY_PROFILES = {
    "acc_tcga":  "acc_tcga_rna_seq_v2_mrna",
    "blca_tcga": "blca_tcga_rna_seq_v2_mrna",
    "brca_tcga": "brca_tcga_rna_seq_v2_mrna",
    "cesc_tcga": "cesc_tcga_rna_seq_v2_mrna",
    "chol_tcga": "chol_tcga_rna_seq_v2_mrna",
    "esca_tcga": "esca_tcga_rna_seq_v2_mrna",
    "gbm_tcga":  "gbm_tcga_rna_seq_v2_mrna",
    "hnsc_tcga": "hnsc_tcga_rna_seq_v2_mrna",
    "kich_tcga": "kich_tcga_rna_seq_v2_mrna",
    "lgg_tcga":  "lgg_tcga_rna_seq_v2_mrna",
    "lihc_tcga": "lihc_tcga_rna_seq_v2_mrna",
    "luad_tcga": "luad_tcga_rna_seq_v2_mrna",
    "lusc_tcga": "lusc_tcga_rna_seq_v2_mrna",
    "ov_tcga":   "ov_tcga_rna_seq_v2_mrna",
    "paad_tcga": "paad_tcga_rna_seq_v2_mrna",
    "prad_tcga": "prad_tcga_rna_seq_v2_mrna",
    "sarc_tcga": "sarc_tcga_rna_seq_v2_mrna",
    "skcm_tcga": "skcm_tcga_rna_seq_v2_mrna",
    "stad_tcga": "stad_tcga_rna_seq_mrna",
    "ucec_tcga": "ucec_tcga_rna_seq_v2_mrna",
    "ucs_tcga":  "ucs_tcga_rna_seq_v2_mrna",
}

def fetch_gene(profile_id, study_id, entrez_id, gene_symbol):
    """Fetch expression for one gene in one study."""
    url = (f"{BASE}/molecular-profiles/{profile_id}/molecular-data"
           f"?sampleListId={study_id}_all&entrezGeneId={entrez_id}")
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200 and r.json():
            return [(e["sampleId"], gene_symbol, e["value"], study_id)
                    for e in r.json() if e.get("value") is not None]
    except Exception:
        pass
    return []

# Fetch all genes for all studies
print(f"Fetching {len(EMT_GENES)} genes x {len(STUDY_PROFILES)} studies...\n")

all_rows = []
for study_id, profile_id in STUDY_PROFILES.items():
    rows_this_study = []
    for gene_symbol, entrez_id in EMT_GENES.items():
        rows = fetch_gene(profile_id, study_id, entrez_id, gene_symbol)
        rows_this_study.extend(rows)
        time.sleep(0.05)  # be gentle with the API

    n_samples = len(set(r[0] for r in rows_this_study))
    n_genes   = len(set(r[1] for r in rows_this_study))
    print(f"  {study_id}: {n_samples} samples, {n_genes} genes")
    all_rows.extend(rows_this_study)

print(f"\nTotal rows: {len(all_rows):,}")

# Build expression matrix
df = pd.DataFrame(all_rows, columns=["sample_id","gene","value","cancer_type"])

expr = df.pivot_table(index=["sample_id","cancer_type"],
                      columns="gene", values="value", aggfunc="mean")
expr = expr.reset_index()

print(f"Expression matrix: {expr.shape[0]} samples x {expr.shape[1]-2} genes")
print(f"Cancer types: {df['cancer_type'].nunique()}")
print(f"Total unique samples: {df['sample_id'].nunique():,}")

out = "data/expression/emt_expression_matrix.tsv"
expr.to_csv(out, sep="\t", index=False)
print(f"\nSaved to: {out}")
