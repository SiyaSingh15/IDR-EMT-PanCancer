import pandas as pd
import os

# ---------------------------------------------------------------
# Build barcode -> cancer type lookup from GDC sample sheet
# ---------------------------------------------------------------
print("Building barcode to cancer type lookup...")

cases = pd.read_csv("data/tcga_sample_cancer_map.tsv", sep="\t")

# Melt all sample ID columns into one long table
id_cols = [c for c in cases.columns if "submitter_sample_ids" in c]
barcode_map = {}

for _, row in cases.iterrows():
    project = row["project.project_id"]  # e.g. TCGA-BRCA
    for col in id_cols:
        barcode = str(row[col])
        if barcode != "nan" and barcode.startswith("TCGA-"):
            # Use first 15 chars as key (matches our sample_id_short format)
            barcode_map[barcode[:15]] = project

print(f"Total barcode mappings: {len(barcode_map):,}")
print(f"Unique cancer types: {len(set(barcode_map.values()))}")

# ---------------------------------------------------------------
# Fix mutation table cancer types
# ---------------------------------------------------------------
muts = pd.read_csv("results/emt_mutations_idr_flagged.tsv", sep="\t")
muts["sample_short"] = muts["Tumor_Sample_Barcode"].str[:15]
muts["cancer_type"] = muts["sample_short"].map(barcode_map).fillna("UNKNOWN")

print(f"\nCorrected cancer types in mutation table:")
print(muts["cancer_type"].value_counts().head(20).to_string())

muts.drop(columns=["sample_short"], inplace=True)
muts.to_csv("results/emt_mutations_idr_flagged.tsv", sep="\t", index=False)
print(f"\nSaved corrected mutations: {len(muts):,} rows")

# ---------------------------------------------------------------
# Fix master dataset cancer types for mutation-flagged samples
# ---------------------------------------------------------------
master = pd.read_csv("results/master_dataset.tsv", sep="\t")
master["cancer_type_from_expr"] = master["cancer_type"].copy()

# Map scored samples using barcode lookup
master["cancer_type_corrected"] = master["sample_id"].str[:15].map(
    barcode_map).fillna(master["cancer_type"])

print(f"\nCorrected cancer types in master dataset:")
print(master["cancer_type_corrected"].value_counts().head(20).to_string())

# Use corrected cancer type
master["cancer_type"] = master["cancer_type_corrected"]
master.drop(columns=["cancer_type_corrected","cancer_type_from_expr"],
            inplace=True)

master.to_csv("results/master_dataset.tsv", sep="\t", index=False)
print(f"\nMaster dataset saved: {len(master):,} rows")

# ---------------------------------------------------------------
# Show IDR mutations per corrected cancer type
# ---------------------------------------------------------------
idr_muts = muts[(muts["in_idr"] == True) &
                (muts["Variant_Classification"] != "Silent")]
print(f"\nIDR mutations by corrected cancer type:")
print(idr_muts["cancer_type"].value_counts().to_string())