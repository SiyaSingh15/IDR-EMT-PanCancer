import pandas as pd
import re
import os

RESULTS_DIR = "results"

mutations = pd.read_csv("results/emt_mutations_all.tsv", sep="\t")
idrs = pd.read_csv("results/idr_boundaries.tsv", sep="\t")

print(f"Mutations loaded: {len(mutations):,}")
print(f"IDR regions loaded: {len(idrs)}")

# Fix cancer type: MC3 barcodes are TCGA-XX-XXXX-01 
# --> Cancer type is encoded in the TCGA project code, not the barcode prefix 
# --> We need a lookup table mapping patient codes to cancer types
# --> MC3 barcode format: TCGA-{2char_tissue_code}-{patient}-{sample}
# --> The 2-letter tissue code IS the cancer type identifier

def extract_cancer_type(barcode):
    """Extract TCGA cancer type from barcode e.g. TCGA-OR-A5J1-01 -> TCGA-ACC"""
   
    parts = str(barcode).split("-")
    if len(parts) >= 2:
        return "TCGA-" + parts[1]
    return "UNKNOWN"

# TCGA 2-letter code to cancer type mapping
TCGA_CODE_MAP = {
    "OR": "TCGA-ACC", "OL": "TCGA-ACC",
    "BT": "TCGA-BLCA", "CU": "TCGA-BLCA", "GC": "TCGA-BLCA",
    "A1": "TCGA-BLCA", "BL": "TCGA-BLCA", "XF": "TCGA-BLCA",
    "A2": "TCGA-BRCA", "AO": "TCGA-BRCA", "AR": "TCGA-BRCA",
    "BH": "TCGA-BRCA", "D8": "TCGA-BRCA", "E2": "TCGA-BRCA",
    "EW": "TCGA-BRCA", "GI": "TCGA-BRCA",
    "C5": "TCGA-CESC", "C6": "TCGA-CESC", "EA": "TCGA-CESC",
    "EX": "TCGA-CESC", "FU": "TCGA-CESC", "HM": "TCGA-CESC",
    "VS": "TCGA-CESC",
    "W5": "TCGA-CHOL", "ZH": "TCGA-CHOL",
    "AA": "TCGA-COAD", "AD": "TCGA-COAD", "AF": "TCGA-COAD",
    "AG": "TCGA-COAD", "AH": "TCGA-COAD", "AY": "TCGA-COAD",
    "CK": "TCGA-COAD", "CM": "TCGA-COAD", "DM": "TCGA-COAD",
    "F4": "TCGA-COAD", "QG": "TCGA-COAD",
    "L5": "TCGA-DLBC", "GS": "TCGA-DLBC",
    "LN": "TCGA-ESCA", "Q9": "TCGA-ESCA", "R6": "TCGA-ESCA",
    "VR": "TCGA-ESCA", "IG": "TCGA-ESCA",
    "02": "TCGA-GBM", "06": "TCGA-GBM", "08": "TCGA-GBM",
    "12": "TCGA-GBM", "14": "TCGA-GBM", "15": "TCGA-GBM",
    "16": "TCGA-GBM", "19": "TCGA-GBM", "32": "TCGA-GBM",
    "4W": "TCGA-GBM",
    "BA": "TCGA-HNSC", "BB": "TCGA-HNSC", "BC": "TCGA-HNSC",
    "BD": "TCGA-HNSC", "CN": "TCGA-HNSC", "CR": "TCGA-HNSC",
    "CV": "TCGA-HNSC", "CX": "TCGA-HNSC", "D6": "TCGA-HNSC",
    "F7": "TCGA-HNSC", "HD": "TCGA-HNSC", "IQ": "TCGA-HNSC",
    "KU": "TCGA-HNSC", "P3": "TCGA-HNSC", "QK": "TCGA-HNSC",
    "T2": "TCGA-HNSC", "TN": "TCGA-HNSC", "UF": "TCGA-HNSC",
    "UL": "TCGA-HNSC",
    "KL": "TCGA-KICH", "KM": "TCGA-KICH", "KN": "TCGA-KICH",
    "B0": "TCGA-KIRC", "B2": "TCGA-KIRC", "B8": "TCGA-KIRC",
    "BP": "TCGA-KIRC", "CJ": "TCGA-KIRC", "CZ": "TCGA-KIRC",
    "DV": "TCGA-KIRC",
    "5P": "TCGA-KIRP", "AK": "TCGA-KIRP", "AL": "TCGA-KIRP",
    "AM": "TCGA-KIRP", "BQ": "TCGA-KIRP", "CZ": "TCGA-KIRP",
    "HE": "TCGA-KIRP", "IA": "TCGA-KIRP", "UW": "TCGA-KIRP",
    "AB": "TCGA-LAML",
    "CS": "TCGA-LGG", "DB": "TCGA-LGG", "DD": "TCGA-LGG",
    "DF": "TCGA-LGG", "DH": "TCGA-LGG", "DU": "TCGA-LGG",
    "DW": "TCGA-LGG", "E1": "TCGA-LGG", "HT": "TCGA-LGG",
    "P5": "TCGA-LGG", "QH": "TCGA-LGG", "R8": "TCGA-LGG",
    "S9": "TCGA-LGG", "TQ": "TCGA-LGG", "VM": "TCGA-LGG",
    "VV": "TCGA-LGG",
    "BC": "TCGA-LIHC", "DD": "TCGA-LIHC", "ES": "TCGA-LIHC",
    "FV": "TCGA-LIHC", "G3": "TCGA-LIHC", "MI": "TCGA-LIHC",
    "NX": "TCGA-LIHC", "RC": "TCGA-LIHC", "RB": "TCGA-LIHC",
    "UB": "TCGA-LIHC", "WQ": "TCGA-LIHC", "YA": "TCGA-LIHC",
    "ZS": "TCGA-LIHC",
    "4B": "TCGA-LUAD", "55": "TCGA-LUAD", "67": "TCGA-LUAD",
    "86": "TCGA-LUAD", "95": "TCGA-LUAD", "CK": "TCGA-LUAD",
    "J2": "TCGA-LUAD", "L4": "TCGA-LUAD", "MP": "TCGA-LUAD",
    "MN": "TCGA-LUAD", "NJ": "TCGA-LUAD", "O1": "TCGA-LUAD",
    "P6": "TCGA-LUAD", "QS": "TCGA-LUAD", "RB": "TCGA-LUAD",
    "S7": "TCGA-LUAD", "UC": "TCGA-LUAD",
    "18": "TCGA-LUSC", "22": "TCGA-LUSC", "33": "TCGA-LUSC",
    "37": "TCGA-LUSC", "43": "TCGA-LUSC", "56": "TCGA-LUSC",
    "58": "TCGA-LUSC", "60": "TCGA-LUSC", "63": "TCGA-LUSC",
    "66": "TCGA-LUSC", "6X": "TCGA-LUSC", "85": "TCGA-LUSC",
    "90": "TCGA-LUSC", "92": "TCGA-LUSC", "94": "TCGA-LUSC",
    "96": "TCGA-LUSC", "98": "TCGA-LUSC", "XC": "TCGA-LUSC",
    "13": "TCGA-OV", "23": "TCGA-OV", "24": "TCGA-OV",
    "25": "TCGA-OV", "29": "TCGA-OV", "30": "TCGA-OV",
    "31": "TCGA-OV", "36": "TCGA-OV", "42": "TCGA-OV",
    "57": "TCGA-OV", "59": "TCGA-OV", "61": "TCGA-OV",
    "IB": "TCGA-PAAD", "F2": "TCGA-PAAD", "FQ": "TCGA-PAAD",
    "H6": "TCGA-PAAD", "HV": "TCGA-PAAD", "HZ": "TCGA-PAAD",
    "L1": "TCGA-PAAD", "LB": "TCGA-PAAD", "M8": "TCGA-PAAD",
    "PZ": "TCGA-PAAD", "Q3": "TCGA-PAAD", "RE": "TCGA-PAAD",
    "S4": "TCGA-PAAD", "US": "TCGA-PAAD", "XD": "TCGA-PAAD",
    "YB": "TCGA-PAAD",
    "CH": "TCGA-PCPG", "E8": "TCGA-PCPG", "QR": "TCGA-PCPG",
    "S2": "TCGA-PCPG", "SR": "TCGA-PCPG", "WM": "TCGA-PCPG",
    "VP": "TCGA-PCPG",
    "CH": "TCGA-PRAD", "EJ": "TCGA-PRAD", "G9": "TCGA-PRAD",
    "GV": "TCGA-PRAD", "H9": "TCGA-PRAD", "HC": "TCGA-PRAD",
    "J4": "TCGA-PRAD", "K7": "TCGA-PRAD", "KC": "TCGA-PRAD",
    "KK": "TCGA-PRAD", "M7": "TCGA-PRAD", "PG": "TCGA-PRAD",
    "V1": "TCGA-PRAD", "XJ": "TCGA-PRAD", "YL": "TCGA-PRAD",
    "ZG": "TCGA-PRAD", "ZJ": "TCGA-PRAD",
    "DC": "TCGA-READ", "DY": "TCGA-READ", "EF": "TCGA-READ",
    "G5": "TCGA-READ",
    "DX": "TCGA-SARC", "MV": "TCGA-SARC", "PB": "TCGA-SARC",
    "PT": "TCGA-SARC", "QC": "TCGA-SARC", "SF": "TCGA-SARC",
    "SY": "TCGA-SARC", "T3": "TCGA-SARC", "UE": "TCGA-SARC",
    "VG": "TCGA-SARC", "WK": "TCGA-SARC", "YN": "TCGA-SARC",
    "DA": "TCGA-SKCM", "EB": "TCGA-SKCM", "ER": "TCGA-SKCM",
    "FD": "TCGA-SKCM", "FE": "TCGA-SKCM", "FJ": "TCGA-SKCM",
    "FP": "TCGA-SKCM", "FR": "TCGA-SKCM", "FS": "TCGA-SKCM",
    "FW": "TCGA-SKCM", "GF": "TCGA-SKCM", "GN": "TCGA-SKCM",
    "HR": "TCGA-SKCM", "HU": "TCGA-SKCM", "RQ": "TCGA-SKCM",
    "WE": "TCGA-SKCM", "YC": "TCGA-SKCM", "YD": "TCGA-SKCM",
    "YG": "TCGA-SKCM", "YN": "TCGA-SKCM", "ZF": "TCGA-SKCM",
    "AY": "TCGA-STAD", "BR": "TCGA-STAD", "CG": "TCGA-STAD",
    "D7": "TCGA-STAD", "F7": "TCGA-STAD", "HB": "TCGA-STAD",
    "KY": "TCGA-STAD", "MQ": "TCGA-STAD", "RD": "TCGA-STAD",
    "SJ": "TCGA-STAD", "VQ": "TCGA-STAD", "XF": "TCGA-STAD",
    "BF": "TCGA-TGCT", "P3": "TCGA-TGCT", "S8": "TCGA-TGCT",
    "2K": "TCGA-THCA", "BJ": "TCGA-THCA", "DE": "TCGA-THCA",
    "DJ": "TCGA-THCA", "DO": "TCGA-THCA", "EL": "TCGA-THCA",
    "EM": "TCGA-THCA", "ET": "TCGA-THCA", "J8": "TCGA-THCA",
    "KS": "TCGA-THCA", "L6": "TCGA-THCA", "QD": "TCGA-THCA",
    "RS": "TCGA-THCA",
    "GF": "TCGA-THYM", "T7": "TCGA-THYM", "XG": "TCGA-THYM",
    "ZT": "TCGA-THYM",
    "AX": "TCGA-UCEC", "AP": "TCGA-UCEC", "BG": "TCGA-UCEC",
    "BS": "TCGA-UCEC", "D1": "TCGA-UCEC", "DI": "TCGA-UCEC",
    "E6": "TCGA-UCEC", "EO": "TCGA-UCEC", "EY": "TCGA-UCEC",
    "FI": "TCGA-UCEC", "GC": "TCGA-UCEC",
    "N1": "TCGA-UCS",
    "RV": "TCGA-UVM", "VD": "TCGA-UVM", "WC": "TCGA-UVM",
    "YZ": "TCGA-UVM",
}

def get_cancer_type(barcode):
    parts = str(barcode).split("-")
    if len(parts) >= 2:
        code = parts[1]
        return TCGA_CODE_MAP.get(code, f"TCGA-{code}")
    return "UNKNOWN"

mutations["cancer_type"] = mutations["Tumor_Sample_Barcode"].apply(get_cancer_type)

print("\nCancer types (mapped):")
print(mutations["cancer_type"].value_counts().head(20).to_string())

# Parse protein position
def parse_position(row):
    pos_str = str(row.get("Protein_position", ""))
    if "/" in pos_str:
        try:
            return int(pos_str.split("/")[0])
        except:
            pass
    hgvs = str(row.get("HGVSp_Short", ""))
    if hgvs.startswith("p."):
        import re
        match = re.search(r'\d+', hgvs)
        if match:
            return int(match.group())
    return None

mutations["protein_pos"] = mutations.apply(parse_position, axis=1)

# Flag IDR mutations
def is_in_idr(gene, pos, idr_df):
    if pd.isna(pos):
        return None
    gene_idrs = idr_df[idr_df["gene"] == gene]
    for _, row in gene_idrs.iterrows():
        if row["idr_start"] <= pos <= row["idr_end"]:
            return True
    return False

def get_idr_region(gene, pos, idr_df):
    if pd.isna(pos):
        return None
    gene_idrs = idr_df[idr_df["gene"] == gene]
    for _, row in gene_idrs.iterrows():
        if row["idr_start"] <= pos <= row["idr_end"]:
            return f"{int(row['idr_start'])}-{int(row['idr_end'])}"
    return None

print("\nFlagging IDR mutations...")
mutations["in_idr"] = mutations.apply(
    lambda r: is_in_idr(r["Hugo_Symbol"], r["protein_pos"], idrs), axis=1)
mutations["idr_region"] = mutations.apply(
    lambda r: get_idr_region(r["Hugo_Symbol"], r["protein_pos"], idrs), axis=1)

coding = mutations[mutations["Variant_Classification"] != "Silent"].copy()
idr_muts = coding[coding["in_idr"] == True]
ord_muts = coding[coding["in_idr"] == False]

print(f"\n=== RESULTS ===")
print(f"Total mutations:           {len(mutations):,}")
print(f"Coding (non-silent):       {len(coding):,}")
print(f"In IDR:                    {len(idr_muts):,} ({100*len(idr_muts)/len(coding):.1f}%)")
print(f"In ordered domain:         {len(ord_muts):,} ({100*len(ord_muts)/len(coding):.1f}%)")

print(f"\nIDR mutations by gene:")
print(idr_muts["Hugo_Symbol"].value_counts().to_string())

print(f"\nIDR mutations by cancer type:")
print(idr_muts["cancer_type"].value_counts().head(15).to_string())

print(f"\nUnique patients with IDR mutations: {idr_muts['Tumor_Sample_Barcode'].nunique():,}")

out = "results/emt_mutations_idr_flagged.tsv"
mutations.to_csv(out, sep="\t", index=False)
print(f"\nSaved to: {out}")
