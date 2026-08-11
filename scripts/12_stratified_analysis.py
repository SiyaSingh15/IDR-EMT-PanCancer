import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import os

RESULTS_DIR = "results"
FIGURES_DIR = "figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

df = pd.read_csv("results/master_dataset.tsv", sep="\t")
print(f"Total samples: {len(df):,}")
print(f"IDR-mutant: {df['has_idr_mutation'].sum()}")

# Stratified analysis: within each cancer type compare IDR-mutant vs IDR-WT EMP scores
print("\n=== WITHIN-CANCER-TYPE ANALYSIS ===\n")

results = []

for cancer in sorted(df["cancer_type"].unique()):
    sub = df[df["cancer_type"] == cancer]
    idr_m = sub[sub["has_idr_mutation"] == True]["emp_score"]
    idr_w = sub[sub["has_idr_mutation"] == False]["emp_score"]

    if len(idr_m) < 3:
        continue

    # Mann-Whitney U test
    stat, pval = stats.mannwhitneyu(idr_m, idr_w, alternative="two-sided")

    # Hybrid enrichment Fisher's exact
    n_hyb_m = (sub[sub["has_idr_mutation"]==True]["emt_state"]=="Hybrid").sum()
    n_oth_m = (sub[sub["has_idr_mutation"]==True]["emt_state"]!="Hybrid").sum()
    n_hyb_w = (sub[sub["has_idr_mutation"]==False]["emt_state"]=="Hybrid").sum()
    n_oth_w = (sub[sub["has_idr_mutation"]==False]["emt_state"]!="Hybrid").sum()
    odds, fp = fisher_exact([[n_hyb_m, n_oth_m],[n_hyb_w, n_oth_w]])

    results.append({
        "cancer_type": cancer,
        "n_total": len(sub),
        "n_idr_mut": len(idr_m),
        "n_idr_wt": len(idr_w),
        "median_emp_idr_mut": idr_m.median(),
        "median_emp_idr_wt": idr_w.median(),
        "delta_emp": idr_m.median() - idr_w.median(),
        "mw_stat": stat,
        "mw_pval": pval,
        "pct_hybrid_idr_mut": 100*n_hyb_m/max(len(idr_m),1),
        "pct_hybrid_idr_wt": 100*n_hyb_w/max(len(idr_w),1),
        "fisher_odds": odds,
        "fisher_pval": fp,
    })

res_df = pd.DataFrame(results)

# FDR correction
if len(res_df) > 0:
    _, mw_fdr, _, _ = multipletests(res_df["mw_pval"], method="fdr_bh")
    _, f_fdr, _, _ = multipletests(res_df["fisher_pval"], method="fdr_bh")
    res_df["mw_fdr"] = mw_fdr
    res_df["fisher_fdr"] = f_fdr

print(res_df[["cancer_type","n_idr_mut","median_emp_idr_mut",
              "median_emp_idr_wt","delta_emp",
              "mw_pval","mw_fdr"]].to_string(index=False))

# Per-gene analysis: ZEB1 and ZEB2 (most mutations)
print("\n=== PER-GENE ANALYSIS ===\n")

gene_results = []
for gene in ["ZEB1","ZEB2","ESRP1","SNAI1","ESRP2"]:
    col = f"idr_mut_{gene}"
    if col not in df.columns:
        continue
    mut = df[df[col]==True]["emp_score"]
    wt  = df[df[col]==False]["emp_score"]
    if len(mut) < 3:
        continue
    stat, pval = stats.mannwhitneyu(mut, wt, alternative="two-sided")
    n_hyb_m = (df[df[col]==True]["emt_state"]=="Hybrid").sum()
    n_oth_m = (df[df[col]==True]["emt_state"]!="Hybrid").sum()
    n_hyb_w = (df[df[col]==False]["emt_state"]=="Hybrid").sum()
    n_oth_w = (df[df[col]==False]["emt_state"]!="Hybrid").sum()
    odds, fp = fisher_exact([[n_hyb_m,n_oth_m],[n_hyb_w,n_oth_w]])
    gene_results.append({
        "gene": gene,
        "n_idr_mut": len(mut),
        "median_emp_mut": mut.median(),
        "median_emp_wt": wt.median(),
        "delta_emp": mut.median()-wt.median(),
        "mw_pval": pval,
        "pct_hybrid_mut": 100*n_hyb_m/max(len(mut),1),
        "pct_hybrid_wt": 100*n_hyb_w/max(len(wt),1),
        "fisher_odds": odds,
        "fisher_pval": fp,
    })

gene_df = pd.DataFrame(gene_results)
if len(gene_df) > 0:
    _, mw_fdr, _, _ = multipletests(gene_df["mw_pval"], method="fdr_bh")
    gene_df["mw_fdr"] = mw_fdr
    print(gene_df[["gene","n_idr_mut","median_emp_mut","median_emp_wt",
                   "delta_emp","mw_pval","mw_fdr",
                   "pct_hybrid_mut","pct_hybrid_wt","fisher_pval"]].to_string(index=False))

# Figure 4: Forest plot of delta EMP scores per cancer type
plot_df = res_df.sort_values("delta_emp").copy()
colors = ["#E74C3C" if p < 0.05 else "#95A5A6" for p in plot_df["mw_pval"]]

fig, ax = plt.subplots(figsize=(8, max(6, len(plot_df)*0.4)))
bars = ax.barh(range(len(plot_df)), plot_df["delta_emp"],
               color=colors, alpha=0.8, height=0.6)
ax.axvline(0, color="black", linewidth=0.8)
ax.set_yticks(range(len(plot_df)))
ax.set_yticklabels([c.replace("_tcga","").replace("TCGA-","")
                    for c in plot_df["cancer_type"]], fontsize=9)
ax.set_xlabel("Delta EMP score (IDR-mutant − IDR-WT median)", fontsize=11)
ax.set_title("IDR mutation effect on EMP score by cancer type\n"
             "(red = p<0.05 nominal, grey = n.s.)", fontsize=11)

# Annotate n
for i, row in enumerate(plot_df.itertuples()):
    ax.text(plot_df["delta_emp"].max()*0.02 if row.delta_emp >= 0
            else plot_df["delta_emp"].min()*0.02,
            i, f"n={row.n_idr_mut}", va="center", fontsize=7,
            color="white" if abs(row.delta_emp) > abs(plot_df["delta_emp"]).max()*0.3 else "black")

plt.tight_layout()
plt.savefig("figures/fig4_forest_delta_emp.png", dpi=150)
plt.close()
print("\nSaved: figures/fig4_forest_delta_emp.png")

# Figure 5: ZEB1 IDR-mutant vs WT violin by cancer type
zeb1_cancers = df[df["idr_mut_ZEB1"]==True]["cancer_type"].value_counts()
zeb1_cancers = zeb1_cancers[zeb1_cancers >= 3].index.tolist()

if zeb1_cancers:
    plot_data = df[df["cancer_type"].isin(zeb1_cancers)].copy()
    plot_data["Group"] = plot_data.apply(
        lambda r: "ZEB1 IDR-mut" if r["idr_mut_ZEB1"] else "ZEB1 IDR-WT", axis=1)
    plot_data["cancer_short"] = plot_data["cancer_type"].str.replace("TCGA-","")

    fig, ax = plt.subplots(figsize=(12, 5))
    sns.boxplot(data=plot_data, x="cancer_short", y="emp_score",
                hue="Group", palette=["#E74C3C","#95A5A6"],
                ax=ax, width=0.6)
    ax.set_xlabel("Cancer type", fontsize=11)
    ax.set_ylabel("EMP score", fontsize=11)
    ax.set_title("ZEB1 IDR mutations vs EMP score (cancer types with n≥3)", fontsize=11)
    ax.legend(title="", fontsize=9)
    plt.tight_layout()
    plt.savefig("figures/fig5_zeb1_idr_emp.png", dpi=150)
    plt.close()
    print("Saved: figures/fig5_zeb1_idr_emp.png")

# Save results
res_df.to_csv("results/stratified_analysis.tsv", sep="\t", index=False)
gene_df.to_csv("results/gene_analysis.tsv", sep="\t", index=False)
print("\nAll results saved.")
