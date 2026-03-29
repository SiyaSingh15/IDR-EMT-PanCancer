import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import os

RESULTS_DIR = "results"
FIGURES_DIR = "figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

# ---------------------------------------------------------------
# Load master dataset
# ---------------------------------------------------------------
df = pd.read_csv("results/master_dataset.tsv", sep="\t")
print(f"Master dataset: {len(df):,} samples")
print(f"IDR-mutant samples: {df['has_idr_mutation'].sum()}")
print(f"Any-EMT-mutant samples: {df['has_any_emt_mutation'].sum()}")

# ---------------------------------------------------------------
# Figure 1: EMP score distribution across cancer types
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(14, 6))
order = df.groupby("cancer_type")["emp_score"].median().sort_values().index
sns.boxplot(data=df, x="cancer_type", y="emp_score",
            order=order, palette="coolwarm", ax=ax)
ax.axhline(0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
ax.set_xticklabels([x.replace("_tcga","").upper()
                    for x in order], rotation=45, ha="right", fontsize=9)
ax.set_xlabel("Cancer type", fontsize=11)
ax.set_ylabel("EMP score (mesenchymal - epithelial)", fontsize=11)
ax.set_title("EMP score distribution across 21 TCGA cancer types", fontsize=13)
plt.tight_layout()
plt.savefig("figures/fig1_emp_scores_by_cancer.png", dpi=150)
plt.close()
print("Saved: figures/fig1_emp_scores_by_cancer.png")

# ---------------------------------------------------------------
# Figure 2: EMT state distribution per cancer type
# ---------------------------------------------------------------
state_counts = df.groupby(["cancer_type","emt_state"]).size().unstack(fill_value=0)
state_pct = state_counts.div(state_counts.sum(axis=1), axis=0) * 100
state_pct = state_pct.sort_values("Mesenchymal", ascending=True)

colors = {"Epithelial": "#4A90D9", "Hybrid": "#F5A623", "Mesenchymal": "#D0021B"}
fig, ax = plt.subplots(figsize=(10, 8))
bottom = np.zeros(len(state_pct))
for state in ["Epithelial", "Hybrid", "Mesenchymal"]:
    if state in state_pct.columns:
        vals = state_pct[state].values
        ax.barh(range(len(state_pct)), vals, left=bottom,
                color=colors[state], label=state, alpha=0.85)
        bottom += vals

ax.set_yticks(range(len(state_pct)))
ax.set_yticklabels([x.replace("_tcga","").upper() for x in state_pct.index],
                    fontsize=9)
ax.set_xlabel("Percentage of samples (%)", fontsize=11)
ax.set_title("E / Hybrid / M state composition by cancer type", fontsize=13)
ax.legend(loc="lower right")
ax.axvline(50, color="black", linestyle="--", linewidth=0.5, alpha=0.4)
plt.tight_layout()
plt.savefig("figures/fig2_emt_state_composition.png", dpi=150)
plt.close()
print("Saved: figures/fig2_emt_state_composition.png")

# ---------------------------------------------------------------
# Figure 3: EMP score — IDR mutant vs WT (all samples)
# ---------------------------------------------------------------
idr_mut  = df[df["has_idr_mutation"] == True]["emp_score"]
idr_wt   = df[df["has_idr_mutation"] == False]["emp_score"]

stat, pval = stats.mannwhitneyu(idr_mut, idr_wt, alternative="two-sided")
print(f"\nMann-Whitney U test (IDR-mut vs IDR-WT):")
print(f"  IDR-mutant:  n={len(idr_mut)}, median={idr_mut.median():.1f}")
print(f"  IDR-WT:      n={len(idr_wt)}, median={idr_wt.median():.1f}")
print(f"  U={stat:.0f}, p={pval:.4f}")

fig, ax = plt.subplots(figsize=(6, 5))
plot_df = pd.DataFrame({
    "EMP score": pd.concat([idr_mut, idr_wt]),
    "Group": ["IDR-mutant"]*len(idr_mut) + ["IDR-WT"]*len(idr_wt)
})
sns.violinplot(data=plot_df, x="Group", y="EMP score",
               palette=["#E74C3C","#95A5A6"], inner="box", ax=ax)
ax.set_title(f"EMP score: IDR-mutant vs IDR-WT\np={pval:.4f} (Mann-Whitney U)",
             fontsize=11)
ax.set_ylabel("EMP score", fontsize=11)
plt.tight_layout()
plt.savefig("figures/fig3_emp_idr_mut_vs_wt.png", dpi=150)
plt.close()
print("Saved: figures/fig3_emp_idr_mut_vs_wt.png")

# ---------------------------------------------------------------
# Figure 4: Hybrid state enrichment — IDR mutant vs WT
# ---------------------------------------------------------------
from scipy.stats import fisher_exact

def hybrid_enrichment(df_sub, label):
    idr_m  = df_sub[df_sub["has_idr_mutation"] == True]
    idr_w  = df_sub[df_sub["has_idr_mutation"] == False]
    if len(idr_m) == 0:
        return None
    n_hybrid_mut = (idr_m["emt_state"] == "Hybrid").sum()
    n_other_mut  = (idr_m["emt_state"] != "Hybrid").sum()
    n_hybrid_wt  = (idr_w["emt_state"] == "Hybrid").sum()
    n_other_wt   = (idr_w["emt_state"] != "Hybrid").sum()
    table = [[n_hybrid_mut, n_other_mut],
             [n_hybrid_wt,  n_other_wt]]
    odds, p = fisher_exact(table, alternative="greater")
    return {
        "label": label,
        "n_idr_mutant": len(idr_m),
        "pct_hybrid_idr_mut": 100*n_hybrid_mut/max(len(idr_m),1),
        "pct_hybrid_idr_wt":  100*n_hybrid_wt/max(len(idr_w),1),
        "odds_ratio": odds,
        "p_value": p
    }

result = hybrid_enrichment(df, "Pan-cancer")
if result:
    print(f"\nFisher's exact test (Hybrid enrichment in IDR-mutant):")
    print(f"  % Hybrid in IDR-mutant: {result['pct_hybrid_idr_mut']:.1f}%")
    print(f"  % Hybrid in IDR-WT:     {result['pct_hybrid_idr_wt']:.1f}%")
    print(f"  Odds ratio: {result['odds_ratio']:.3f}, p={result['p_value']:.4f}")

# ---------------------------------------------------------------
# Save summary stats
# ---------------------------------------------------------------
summary = {
    "total_samples": len(df),
    "idr_mutant_samples": int(df["has_idr_mutation"].sum()),
    "any_emt_mutant_samples": int(df["has_any_emt_mutation"].sum()),
    "mannwhitney_U": float(stat),
    "mannwhitney_p": float(pval),
    "idr_mut_median_emp": float(idr_mut.median()) if len(idr_mut) > 0 else None,
    "idr_wt_median_emp": float(idr_wt.median()),
}
import json
with open("results/stats_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print(f"\nStats summary saved to: results/stats_summary.json")
print("\nWeek 2 complete! All figures saved to figures/")