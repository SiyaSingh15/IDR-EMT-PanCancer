import requests
import json
import pandas as pd
import os

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

EMT_PROTEINS = {
    "SNAI1":  "O95863",
    "SNAI2":  "O43623",
    "ZEB1":   "P37275",
    "ZEB2":   "O60315",
    "TWIST1": "Q15672",
    "TWIST2": "Q8WVJ9",
    "ESRP1":  "Q9NWF9",
    "ESRP2":  "Q9UBB5",
    "CDH1":   "P12830",
    "VIM":    "P08670",
}

# We use curated-disorder-merge as PRIMARY source (experimentally supported)
# and prediction-disorder-mobidb_lite as SECONDARY (computational consensus)
PRIMARY_KEY   = "curated-disorder-merge"
SECONDARY_KEY = "prediction-disorder-mobidb_lite"

def fetch_mobidb(uniprot_id):
    url = f"https://mobidb.org/{uniprot_id}?format=json"
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200:
            return r.json()
    except Exception as e:
        print(f"Error: {e}")
    return None

def extract_idrs(data, gene):
    """
    Extract IDR regions using curated data first, 
    falling back to mobidb_lite prediction.
    Returns list of (start, end) tuples and the source used.
    """
    protein_length = data.get("length", "?")

    # Try curated first
    for key in [PRIMARY_KEY, SECONDARY_KEY]:
        ann = data.get(key)
        if ann and "regions" in ann:
            regions = [(int(r[0]), int(r[1])) for r in ann["regions"]]
            return regions, key, protein_length

    return [], "none", protein_length

# ---------------------------------------------------------------
# Load already-fetched raw data (no need to re-query API)
# ---------------------------------------------------------------
print("Loading MobiDB data from saved JSON...\n")

with open("results/mobidb_raw.json") as f:
    raw = json.load(f)

all_rows = []

for gene, uid in EMT_PROTEINS.items():
    data = raw.get(gene)
    if not data:
        print(f"  {gene}: no data")
        continue

    idrs, source, prot_len = extract_idrs(data, gene)
    total_idr_res = sum(e - s + 1 for s, e in idrs)

    print(f"  {gene} (len={prot_len}): {len(idrs)} IDR regions, "
          f"{total_idr_res} disordered residues [source: {source}]")

    for start, end in idrs:
        all_rows.append({
            "gene":        gene,
            "uniprot_id":  uid,
            "idr_start":   start,
            "idr_end":     end,
            "idr_length":  end - start + 1,
            "protein_length": prot_len,
            "source":      source
        })

# ---------------------------------------------------------------
# Save and summarise
# ---------------------------------------------------------------
idr_df = pd.DataFrame(all_rows)
idr_df.to_csv("results/idr_boundaries.tsv", sep="\t", index=False)

print(f"\n=== IDR SUMMARY ===")
print(f"Total IDR regions across all genes: {len(idr_df)}")
print()
print(idr_df.groupby("gene").agg(
    protein_length=("protein_length", "first"),
    n_idr_regions=("idr_start", "count"),
    total_idr_residues=("idr_length", "sum"),
    idr_start=("idr_start", list),
    idr_end=("idr_end", list)
).to_string())

print(f"\nSaved to: results/idr_boundaries.tsv")