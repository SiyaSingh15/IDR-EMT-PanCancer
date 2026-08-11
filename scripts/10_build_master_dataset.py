import pandas as pd
import os

RESULTS_DIR = "results"
 
# Load both datasets
print("Loading datasets...")

# IDR-flagged mutations 
muts = pd.read_csv("results/emt_mutations_idr_flagged.tsv", sep="\t")
print(f"Mutations: {len(muts)} rows, {muts['Tumor_Sample_Barcode'].nunique()} unique samples")

# EMT scores 
scores = pd.read_csv("results/emt_scores.tsv", sep="\t")
print(f"EMT scores: {len(scores)} samples")

# TCGA mutation barcodes look like: TCGA-XX-XXXX-01A-...
# cBioPortal sample IDs look like:  TCGA-XX-XXXX-01
# We need to match on the first 15 characters (patient + sample type)
print("\nHarmonising sample IDs...")

muts["sample_id_short"] = muts["Tumor_Sample_Barcode"].str[:15]
scores["sample_id_short"] = scores["sample_id"].str[:15]

print(f"Mutation sample IDs (example): {muts['sample_id_short'].iloc[0]}")
print(f"Score sample IDs (example):    {scores['sample_id_short'].iloc[0]}")

# For each scored sample --> flag whether it has an IDR mutation in any of our 10 EMT genes

# Get IDR-mutant samples 
idr_mut_samples = set(
    muts[
        (muts["in_idr"] == True) &
        (muts["Variant_Classification"] != "Silent")
    ]["sample_id_short"]
)

# Get any EMT-mutant samples 
any_mut_samples = set(
    muts[
        muts["Variant_Classification"] != "Silent"
    ]["sample_id_short"]
)

print(f"\nSamples with IDR mutations: {len(idr_mut_samples)}")
print(f"Samples with any EMT coding mutation: {len(any_mut_samples)}")

# Flag each scored sample
scores["has_idr_mutation"]  = scores["sample_id_short"].isin(idr_mut_samples)
scores["has_any_emt_mutation"] = scores["sample_id_short"].isin(any_mut_samples)

# Also add per-gene IDR mutation flags
for gene in ["SNAI1","SNAI2","ZEB1","ZEB2","TWIST1","TWIST2",
             "ESRP1","ESRP2","CDH1","VIM"]:
    gene_idr_samples = set(
        muts[
            (muts["Hugo_Symbol"] == gene) &
            (muts["in_idr"] == True) &
            (muts["Variant_Classification"] != "Silent")
        ]["sample_id_short"]
    )
    scores[f"idr_mut_{gene}"] = scores["sample_id_short"].isin(gene_idr_samples)

# Summary statistics
print(f"\n=== MASTER DATASET - SUMMARY ===")
print(f"Total samples: {len(scores):,}")
print(f"Samples with IDR mutation: {scores['has_idr_mutation'].sum()}")
print(f"Samples with any EMT mutation: {scores['has_any_emt_mutation'].sum()}")

print(f"\nEMT state in IDR-mutant vs IDR-WT samples:")
print(scores.groupby("has_idr_mutation")["emt_state"].value_counts().unstack(fill_value=0))

print(f"\nEMT state in any-EMT-mutant vs WT samples:")
print(scores.groupby("has_any_emt_mutation")["emt_state"].value_counts().unstack(fill_value=0))

print(f"\nPer-gene IDR mutation counts:")
for gene in ["SNAI1","SNAI2","ZEB1","ZEB2","TWIST1","TWIST2",
             "ESRP1","ESRP2","CDH1","VIM"]:
    n = scores[f"idr_mut_{gene}"].sum()
    if n > 0:
        print(f"  {gene}: {n} samples with IDR mutation")

# Save master dataset
out = "results/master_dataset.tsv"
scores.to_csv(out, sep="\t", index=False)
print(f"\nMaster dataset saved to: {out}")
print(f"Columns: {list(scores.columns)}")
