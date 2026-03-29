import pandas as pd
import gzip
import os

DATA_DIR = "data/mafs"

# ---------------------------------------------------------------
# Step 1: Check all files are present and show sizes
# ---------------------------------------------------------------
print("=== Downloaded MAF files ===\n")
files = sorted(os.listdir(DATA_DIR))
total_size = 0
for f in files:
    path = os.path.join(DATA_DIR, f)
    size = os.path.getsize(path) / 1e6
    total_size += size
    print(f"  {f}: {size:.1f} MB")
print(f"\nTotal: {total_size:.1f} MB across {len(files)} files")

# ---------------------------------------------------------------
# Step 2: Peek inside one MAF file (TCGA-BRCA as example)
# ---------------------------------------------------------------
print("\n=== Peeking inside TCGA-BRCA MAF ===\n")

brca_path = os.path.join(DATA_DIR, "TCGA-BRCA.maf.gz")

# MAF files have comment lines starting with # at the top
# We skip those and read the actual table
df = pd.read_csv(brca_path, sep="\t", comment="#",
                 low_memory=False, nrows=5)

print(f"Columns in MAF file: {len(df.columns)}")
print(f"\nFirst 5 column names:")
for col in df.columns[:10]:
    print(f"  {col}")

print(f"\nShape of first 5 rows: {df.shape}")

# ---------------------------------------------------------------
# Step 3: Check our EMT genes are present
# ---------------------------------------------------------------
print("\n=== Checking EMT gene presence in TCGA-BRCA ===\n")

# Read full BRCA MAF
df_full = pd.read_csv(brca_path, sep="\t", comment="#", low_memory=False)
print(f"Total mutations in BRCA: {len(df_full):,}")

EMT_GENES = ["SNAI1", "SNAI2", "ZEB1", "ZEB2",
             "TWIST1", "TWIST2", "ESRP1", "ESRP2",
             "CDH1", "VIM"]

print(f"\nMutations in our EMT genes:")
for gene in EMT_GENES:
    count = (df_full["Hugo_Symbol"] == gene).sum()
    print(f"  {gene}: {count} mutations")