import pandas as pd
import os

DATA_DIR = "data/mafs"
RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

EMT_GENES = ["SNAI1", "SNAI2", "ZEB1", "ZEB2",
             "TWIST1", "TWIST2", "ESRP1", "ESRP2",
             "CDH1", "VIM"]

# Columns we need from MC3
COLS = ["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
        "Variant_Classification", "Variant_Type", "Reference_Allele",
        "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "HGVSp_Short",
        "Protein_position", "SWISSPROT", "FILTER", "IMPACT"]

mc3_path = "data/mafs/mc3.v0.2.8.PUBLIC.maf.gz"
print(f"Reading MC3: {os.path.getsize(mc3_path)/1e6:.0f} MB\n")

# Read only the columns we need — much faster and less memory
df = pd.read_csv(mc3_path, sep="\t", low_memory=False,
                 usecols=lambda c: c in COLS)

print(f"Total mutations: {len(df):,}")
print(f"Unique samples:  {df['Tumor_Sample_Barcode'].nunique():,}")

# MC3 recommended filter: keep PASS mutations only
# (removes low-confidence calls flagged by quality filters)
if "FILTER" in df.columns:
    before = len(df)
    df = df[df["FILTER"] == "PASS"].copy()
    print(f"After PASS filter: {len(df):,} ({before - len(df):,} removed)")

# Filter for our EMT genes
df_emt = df[df["Hugo_Symbol"].isin(EMT_GENES)].copy()

# Extract cancer type from barcode (TCGA-XX-...)
df_emt["cancer_type"] = df_emt["Tumor_Sample_Barcode"].str.extract(r'(TCGA-[A-Z]+)')

print(f"\nEMT gene mutations: {len(df_emt):,}")
print(f"Unique patients:    {df_emt['Tumor_Sample_Barcode'].nunique():,}")

print("\nMutations per EMT gene:")
print(df_emt["Hugo_Symbol"].value_counts().to_string())

print("\nVariant classifications:")
print(df_emt["Variant_Classification"].value_counts().to_string())

print("\nCancer types represented:")
print(df_emt["cancer_type"].value_counts().to_string())

out = os.path.join(RESULTS_DIR, "emt_mutations_all.tsv")
df_emt.to_csv(out, sep="\t", index=False)
print(f"\nSaved to: {out}")