import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from scipy.stats import fisher_exact
import os

FIGURES_DIR = "figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

master  = pd.read_csv("results/master_dataset.tsv",    sep="\t")
surv    = pd.read_csv("results/survival_dataset.tsv",  sep="\t")
idr_df  = pd.read_csv("results/idr_boundaries.tsv",    sep="\t")
muts    = pd.read_csv("results/emt_mutations_idr_flagged.tsv", sep="\t")
strat   = pd.read_csv("results/stratified_analysis.tsv", sep="\t")
genes   = pd.read_csv("results/gene_analysis.tsv",     sep="\t")

# Figure A: IDR coverage lollipop — all 10 genes
print("Generating Figure A: IDR lollipop plots...")

fig, axes = plt.subplots(5, 2, figsize=(14, 18))
axes = axes.flatten()

EMT_PROTEINS = {
    "SNAI1": "O95863", "SNAI2": "O43623", "ZEB1": "P37275",
    "ZEB2": "O60315",  "TWIST1": "Q15672", "TWIST2": "Q8WVJ9",
    "ESRP1": "Q9NWF9", "ESRP2": "Q9UBB5", "CDH1": "P12830",
    "VIM": "P08670"
}

PROTEIN_LENGTHS = {
    "SNAI1": 264, "SNAI2": 268, "ZEB1": 1124, "ZEB2": 1214,
    "TWIST1": 202, "TWIST2": 160, "ESRP1": 866, "ESRP2": 411,
    "CDH1": 882, "VIM": 466
}

coding = muts[muts["Variant_Classification"] != "Silent"].copy()
idr_coding = coding[coding["in_idr"] == True]

for idx, gene in enumerate(EMT_PROTEINS.keys()):
    ax = axes[idx]
    plen = PROTEIN_LENGTHS[gene]
    gene_idrs = idr_df[idr_df["gene"] == gene]
    gene_muts = coding[coding["Hugo_Symbol"] == gene]
    gene_idr_muts = idr_coding[idr_coding["Hugo_Symbol"] == gene]

    # Draw protein backbone
    ax.hlines(0, 0, plen, color="#CCCCCC", linewidth=4, zorder=1)

    # Draw IDR regions
    for _, idr_row in gene_idrs.iterrows():
        ax.hlines(0, idr_row["idr_start"], idr_row["idr_end"],
                  color="#3498DB", linewidth=8, alpha=0.6, zorder=2)

    # Plot ordered domain mutations
    ord_muts = gene_muts[gene_muts["in_idr"] != True]
    if len(ord_muts) > 0 and "protein_pos" in ord_muts.columns:
        ord_pos = ord_muts["protein_pos"].dropna()
        ax.vlines(ord_pos, 0, 0.6, color="#AAAAAA", linewidth=0.8, alpha=0.5)
        ax.scatter(ord_pos, [0.6]*len(ord_pos),
                   color="#AAAAAA", s=20, zorder=3, alpha=0.5)

    # Plot IDR mutations
    if len(gene_idr_muts) > 0 and "protein_pos" in gene_idr_muts.columns:
        idr_pos = gene_idr_muts["protein_pos"].dropna()
        nonsense = gene_idr_muts[
            gene_idr_muts["Variant_Classification"]=="Nonsense_Mutation"
        ]["protein_pos"].dropna()
        ax.vlines(idr_pos, 0, 1.0, color="#E74C3C", linewidth=1.2, zorder=4)
        ax.scatter(idr_pos, [1.0]*len(idr_pos),
                   color="#E74C3C", s=35, zorder=5)
        if len(nonsense) > 0:
            ax.scatter(nonsense, [1.0]*len(nonsense),
                       color="#8B0000", s=60, marker="X", zorder=6,
                       label="Nonsense")

    n_idr = len(gene_idr_muts)
    n_total = len(gene_muts)
    ax.set_xlim(-plen*0.02, plen*1.05)
    ax.set_ylim(-0.3, 1.4)
    ax.set_xlabel("Protein position (aa)", fontsize=8)
    ax.set_title(f"{gene}  (n={n_idr} IDR / {n_total} total coding)", fontsize=10)
    ax.set_yticks([])
    ax.spines[["top","right","left"]].set_visible(False)

    # Legend patches
    idr_patch = mpatches.Patch(color="#3498DB", alpha=0.6, label="IDR region")
    mut_patch  = mpatches.Patch(color="#E74C3C", label="IDR mutation")
    ord_patch  = mpatches.Patch(color="#AAAAAA", alpha=0.5, label="Ordered mutation")
    ax.legend(handles=[idr_patch, mut_patch, ord_patch],
              fontsize=6, loc="upper right")

plt.suptitle("Somatic mutations mapped onto IDR boundaries — 10 EMT regulators",
             fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig("figures/figA_lollipop_all_genes.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved: figures/figA_lollipop_all_genes.png")

# Figure B: Stacked summary - IDR vs ordered mutations per gene
print("Generating Figure B: mutation summary...")

gene_summary = []
for gene in EMT_PROTEINS.keys():
    g_muts = coding[coding["Hugo_Symbol"] == gene]
    n_idr  = (g_muts["in_idr"] == True).sum()
    n_ord  = (g_muts["in_idr"] == False).sum()
    n_unk  = g_muts["in_idr"].isna().sum()
    plen   = PROTEIN_LENGTHS[gene]
    idr_len = idr_df[idr_df["gene"]==gene]["idr_length"].sum()
    pct_idr_seq = 100 * idr_len / plen
    pct_idr_mut = 100 * n_idr / max(n_idr+n_ord, 1)
    gene_summary.append({
        "gene": gene, "n_idr_mut": n_idr, "n_ord_mut": n_ord,
        "pct_idr_seq": pct_idr_seq, "pct_idr_mut": pct_idr_mut
    })

gs = pd.DataFrame(gene_summary).sort_values("n_idr_mut", ascending=True)

fig, ax = plt.subplots(figsize=(9, 6))
y = range(len(gs))
ax.barh(y, gs["n_ord_mut"], color="#AAAAAA", alpha=0.7, label="Ordered domain")
ax.barh(y, gs["n_idr_mut"], left=gs["n_ord_mut"],
        color="#E74C3C", alpha=0.85, label="IDR")

ax.set_yticks(y)
ax.set_yticklabels(gs["gene"], fontsize=11)
ax.set_xlabel("Number of coding mutations", fontsize=11)
ax.set_title("Coding mutations in IDR vs ordered domains — 10 EMT regulators\n(MC3 pan-cancer, n=10,295 samples)", fontsize=11)
ax.legend(fontsize=10)

# Annotate % IDR
for i, row in enumerate(gs.itertuples()):
    total = row.n_idr_mut + row.n_ord_mut
    if total > 0:
        ax.text(total + 2, i,
                f"{row.pct_idr_mut:.0f}% IDR ({row.pct_idr_seq:.0f}% seq)",
                va="center", fontsize=8, color="#333333")

ax.set_xlim(0, gs[["n_idr_mut","n_ord_mut"]].sum(axis=1).max() * 1.35)
plt.tight_layout()
plt.savefig("figures/figB_mutation_summary.png", dpi=150)
plt.close()
print("Saved: figures/figB_mutation_summary.png")

# Figure C: ESRP1 IDR mutation - Hybrid state depletion
print("Generating Figure C: ESRP1 hybrid state...")

esrp1_mut = master[master["idr_mut_ESRP1"]==True]
esrp1_wt  = master[master["idr_mut_ESRP1"]==False]

states = ["Epithelial","Hybrid","Mesenchymal"]
mut_pcts = [100*(esrp1_mut["emt_state"]==s).mean() for s in states]
wt_pcts  = [100*(esrp1_wt["emt_state"]==s).mean()  for s in states]

x = np.arange(len(states))
w = 0.35
colors_states = ["#3498DB","#F39C12","#E74C3C"]

fig, ax = plt.subplots(figsize=(7, 5))
bars1 = ax.bar(x-w/2, wt_pcts,  w, label=f"ESRP1 IDR-WT (n={len(esrp1_wt):,})",
               color=colors_states, alpha=0.5)
bars2 = ax.bar(x+w/2, mut_pcts, w, label=f"ESRP1 IDR-mut (n={len(esrp1_mut)})",
               color=colors_states, alpha=0.95)

ax.set_xticks(x)
ax.set_xticklabels(states, fontsize=11)
ax.set_ylabel("% of samples", fontsize=11)
ax.set_title("ESRP1 IDR mutations: depletion of Hybrid E/M state\n"
             f"(Fisher's p=0.096, Hybrid: 28% mut vs 50% WT)", fontsize=11)

# Add value labels
for bar in bars1:
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5,
            f"{bar.get_height():.1f}%", ha="center", va="bottom", fontsize=9)
for bar in bars2:
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5,
            f"{bar.get_height():.1f}%", ha="center", va="bottom",
            fontsize=9, fontweight="bold")

handles = [mpatches.Patch(color=c, label=s)
           for c, s in zip(colors_states, states)]
handles += [mpatches.Patch(color="white",alpha=0,
                            label=f"IDR-WT (lighter) vs IDR-mut (darker)")]
ax.legend(handles=handles[:3], title="EMT state", fontsize=9)
plt.tight_layout()
plt.savefig("figures/figC_esrp1_hybrid_depletion.png", dpi=150)
plt.close()
print("Saved: figures/figC_esrp1_hybrid_depletion.png")

# Figure D: Pan-cancer heatmap - IDR mutation frequency
print("Generating Figure D: pan-cancer heatmap...")

# Count IDR mutations per gene per cancer type
idr_muts = muts[(muts["in_idr"]==True) &
                (muts["Variant_Classification"]!="Silent")]

heatmap_data = idr_muts.groupby(
    ["cancer_type","Hugo_Symbol"]).size().unstack(fill_value=0)

# Keep cancer types with at least 1 IDR mutation
heatmap_data = heatmap_data[heatmap_data.sum(axis=1) > 0]

# Sort by total IDR mutations
heatmap_data = heatmap_data.loc[
    heatmap_data.sum(axis=1).sort_values(ascending=False).index]

fig, ax = plt.subplots(figsize=(12, 8))
sns.heatmap(heatmap_data, annot=True, fmt="d", cmap="YlOrRd",
            linewidths=0.5, ax=ax, cbar_kws={"label":"IDR mutation count"})
ax.set_xlabel("EMT gene", fontsize=11)
ax.set_ylabel("Cancer type", fontsize=11)
ax.set_title("Pan-cancer IDR mutation frequency across 10 EMT regulators\n"
             "(MC3 dataset, PASS-filtered, n=10,295 TCGA samples)", fontsize=12)
plt.tight_layout()
plt.savefig("figures/figD_pancancer_heatmap.png", dpi=150)
plt.close()
print("Saved: figures/figD_pancancer_heatmap.png")

# Print final summary stats 
print("\n" + "="*60)
print("FINAL RESULTS SUMMARY FOR MANUSCRIPT")
print("="*60)

total_muts     = len(coding)
total_idr_muts = (coding["in_idr"]==True).sum()
pct_idr        = 100*total_idr_muts/total_muts

print(f"\nDataset: MC3 pan-cancer (PASS-filtered)")
print(f"Total TCGA samples:          10,295")
print(f"Total coding EMT mutations:  {total_muts:,}")
print(f"IDR-resident mutations:      {total_idr_muts:,} ({pct_idr:.1f}%)")
print(f"Unique patients (IDR-mut):   {muts[muts['in_idr']==True]['Tumor_Sample_Barcode'].nunique():,}")

print(f"\nEMT scoring (cBioPortal RNA-seq, n={len(master):,} samples):")
print(f"  Epithelial:   {(master['emt_state']=='Epithelial').sum():,} ({100*(master['emt_state']=='Epithelial').mean():.1f}%)")
print(f"  Hybrid E/M:   {(master['emt_state']=='Hybrid').sum():,} ({100*(master['emt_state']=='Hybrid').mean():.1f}%)")
print(f"  Mesenchymal:  {(master['emt_state']=='Mesenchymal').sum():,} ({100*(master['emt_state']=='Mesenchymal').mean():.1f}%)")

print(f"\nKey statistical results:")
print(f"  Mann-Whitney (pan-cancer IDR-mut vs WT): p=0.848 (n.s.)")
print(f"  ESRP1 Hybrid depletion Fisher's:         p=0.096")
print(f"  BRCA IDR-mut delta EMP:                  +6,308 (p=0.088)")
print(f"  STAD IDR-mut Mann-Whitney:               p=0.015 (nominal)")
print(f"  Cox HR (IDR mutation):                   0.82 (95%CI 0.64-1.04), p=0.106")
print(f"  Cox HR (Mesenchymal state):              1.20 (95%CI 1.09-1.32), p<0.005")

print("\nAll figures saved to figures/")
