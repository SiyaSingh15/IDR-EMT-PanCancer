import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

FIGURES_DIR = "figures"
RESULTS_DIR = "results"
os.makedirs(FIGURES_DIR, exist_ok=True)

# Build survival dataset
print("Building survival dataset...")

# Patient clinical — pivot to wide
pat_raw = pd.read_csv("data/tcga_patient_clinical.tsv", sep="\t")
pat_wide = pat_raw.pivot_table(
    index=["patientId","studyId"],
    columns="clinicalAttributeId",
    values="value",
    aggfunc="first"
).reset_index()

# Extract OS columns
pat_wide["os_months"] = pd.to_numeric(
    pat_wide.get("OS_MONTHS", pd.Series(dtype=float)), errors="coerce")
pat_wide["os_status"] = pat_wide.get("OS_STATUS","").astype(str).str.upper()
pat_wide["os_event"]  = pat_wide["os_status"].apply(
    lambda x: 1 if "DECEASED" in x or "1:" in x or "DEAD" in x else 0)
pat_wide["age"] = pd.to_numeric(
    pat_wide.get("AGE", pd.Series(dtype=float)), errors="coerce")
pat_wide["cancer_type"] = pat_wide["studyId"].str.upper().str.replace("_TCGA","").apply(
    lambda x: f"TCGA-{x}")

# Keep only patients with valid OS
surv = pat_wide[
    pat_wide["os_months"].notna() &
    (pat_wide["os_months"] > 0)
][["patientId","cancer_type","os_months","os_event","age"]].copy()

print(f"Patients with OS data: {len(surv):,}")
print(f"Events (deaths): {surv['os_event'].sum():,} ({100*surv['os_event'].mean():.1f}%)")

# Load master dataset with IDR mutation flags
master = pd.read_csv("results/master_dataset.tsv", sep="\t")

# Extract patient ID from sample ID (first 12 chars)
master["patientId"] = master["sample_id"].str[:12]

# Merge with survival
merged = surv.merge(
    master[["patientId","emt_state","emp_score","has_idr_mutation",
            "has_any_emt_mutation","idr_mut_ZEB1","idr_mut_ZEB2",
            "idr_mut_ESRP1","idr_mut_SNAI1","cancer_type"]],
    on="patientId", how="inner", suffixes=("","_expr")
)

print(f"\nMerged survival + expression: {len(merged):,} patients")
print(f"IDR-mutant patients: {merged['has_idr_mutation'].sum()}")
print(f"Cancer types: {merged['cancer_type'].nunique()}")

# Figure 6: Pan-cancer KM — IDR-mutant vs IDR-WT
kmf = KaplanMeierFitter()
fig, ax = plt.subplots(figsize=(8, 5))

for label, group, color in [
    ("IDR-mutant", merged[merged["has_idr_mutation"]==True], "#E74C3C"),
    ("IDR-WT",     merged[merged["has_idr_mutation"]==False], "#2980B9")
]:
    if len(group) < 5: continue
    kmf.fit(group["os_months"], event_observed=group["os_event"], label=label)
    kmf.plot_survival_function(ax=ax, color=color, ci_show=True)

# Log-rank test
idr_m = merged[merged["has_idr_mutation"]==True]
idr_w = merged[merged["has_idr_mutation"]==False]
if len(idr_m) >= 5:
    result = logrank_test(
        idr_m["os_months"], idr_w["os_months"],
        event_observed_A=idr_m["os_event"],
        event_observed_B=idr_w["os_event"]
    )
    ax.set_title(f"Overall survival: IDR-mutant vs IDR-WT\n"
                 f"(n_mut={len(idr_m)}, n_wt={len(idr_w)}, "
                 f"log-rank p={result.p_value:.4f})", fontsize=11)
else:
    ax.set_title("Overall survival: IDR-mutant vs IDR-WT", fontsize=11)

ax.set_xlabel("Time (months)", fontsize=11)
ax.set_ylabel("Survival probability", fontsize=11)
ax.set_ylim(0, 1.05)
plt.tight_layout()
plt.savefig("figures/fig6_km_pancancer_idr.png", dpi=150)
plt.close()
print("\nSaved: figures/fig6_km_pancancer_idr.png")

# Figure 7: KM by EMT state (pan-cancer)
fig, ax = plt.subplots(figsize=(8, 5))
colors = {"Epithelial":"#3498DB","Hybrid":"#F39C12","Mesenchymal":"#E74C3C"}
lr_results = {}

for state in ["Epithelial","Hybrid","Mesenchymal"]:
    group = merged[merged["emt_state"]==state]
    if len(group) < 10: continue
    kmf.fit(group["os_months"], event_observed=group["os_event"], label=state)
    kmf.plot_survival_function(ax=ax, color=colors[state], ci_show=False)

ax.set_title("Overall survival by EMT state (pan-cancer)", fontsize=11)
ax.set_xlabel("Time (months)", fontsize=11)
ax.set_ylabel("Survival probability", fontsize=11)
ax.set_ylim(0, 1.05)
plt.tight_layout()
plt.savefig("figures/fig7_km_emt_state.png", dpi=150)
plt.close()
print("Saved: figures/fig7_km_emt_state.png")

# Cox PH model: pan-cancer
print("\n=== COX PROPORTIONAL HAZARDS MODEL ===\n")

cox_data = merged[
    merged["os_months"].notna() &
    merged["age"].notna()
][["os_months","os_event","has_idr_mutation",
   "age","emt_state","cancer_type"]].copy()

# Encode categorical variables
cox_data["emt_mesenchymal"] = (cox_data["emt_state"] == "Mesenchymal").astype(int)
cox_data["emt_hybrid"]      = (cox_data["emt_state"] == "Hybrid").astype(int)
cox_data["idr_mut"]         = cox_data["has_idr_mutation"].astype(int)
cox_data["age_scaled"]      = (cox_data["age"] - cox_data["age"].mean()) / cox_data["age"].std()

# Add cancer type dummies (use BRCA as reference - largest group)
top_cancers = cox_data["cancer_type"].value_counts().head(8).index.tolist()
for ct in top_cancers[1:]:  # skip first as reference
    cox_data[f"ct_{ct}"] = (cox_data["cancer_type"] == ct).astype(int)

ct_cols = [f"ct_{ct}" for ct in top_cancers[1:]]
cox_cols = ["os_months","os_event","idr_mut","age_scaled",
            "emt_mesenchymal","emt_hybrid"] + ct_cols

cox_fit_data = cox_data[cox_cols].dropna()
print(f"Cox model: {len(cox_fit_data):,} patients\n")

try:
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(cox_fit_data, duration_col="os_months", event_col="os_event")
    cph.print_summary()

    # Save forest plot
    fig, ax = plt.subplots(figsize=(8, 6))
    cph.plot(ax=ax)
    ax.set_title("Cox PH model: hazard ratios\n(IDR mutation + EMT state + age + cancer type)",
                 fontsize=11)
    ax.axvline(0, color="black", linewidth=0.8)
    plt.tight_layout()
    plt.savefig("figures/fig8_cox_forest.png", dpi=150)
    plt.close()
    print("\nSaved: figures/fig8_cox_forest.png")

    # Extract IDR mutation hazard ratio
    summary = cph.summary
    if "idr_mut" in summary.index:
        hr = np.exp(summary.loc["idr_mut", "coef"])
        ci_low = np.exp(summary.loc["idr_mut", "coef lower 95%"])
        ci_high = np.exp(summary.loc["idr_mut", "coef upper 95%"])
        pval = summary.loc["idr_mut", "p"]
        print(f"\n>>> IDR mutation hazard ratio: {hr:.3f} "
              f"(95% CI: {ci_low:.3f}-{ci_high:.3f}), p={pval:.4f}")

except Exception as e:
    print(f"Cox model error: {e}")

# Per-cancer KM for top hits: UCEC, LUAD, SKCM
print("\nGenerating per-cancer KM plots...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
top_cancers_km = ["TCGA-UCEC", "TCGA-LUAD", "TCGA-SKCM"]

for ax, cancer in zip(axes, top_cancers_km):
    sub = merged[merged["cancer_type"] == cancer]
    idr_m = sub[sub["has_idr_mutation"]==True]
    idr_w = sub[sub["has_idr_mutation"]==False]

    if len(idr_m) < 3:
        ax.text(0.5, 0.5, f"{cancer}\nn_mut too small",
                ha="center", va="center", transform=ax.transAxes)
        continue

    kmf.fit(idr_m["os_months"], event_observed=idr_m["os_event"],
            label=f"IDR-mut (n={len(idr_m)})")
    kmf.plot_survival_function(ax=ax, color="#E74C3C", ci_show=True)

    kmf.fit(idr_w["os_months"], event_observed=idr_w["os_event"],
            label=f"IDR-WT (n={len(idr_w)})")
    kmf.plot_survival_function(ax=ax, color="#2980B9", ci_show=True)

    result = logrank_test(
        idr_m["os_months"], idr_w["os_months"],
        event_observed_A=idr_m["os_event"],
        event_observed_B=idr_w["os_event"]
    )
    ax.set_title(f"{cancer}\nlog-rank p={result.p_value:.3f}", fontsize=10)
    ax.set_xlabel("Months", fontsize=9)
    ax.set_ylabel("Survival", fontsize=9)
    ax.set_ylim(0, 1.05)

plt.suptitle("Overall survival: IDR-mutant vs IDR-WT by cancer type", fontsize=12)
plt.tight_layout()
plt.savefig("figures/fig9_km_per_cancer.png", dpi=150)
plt.close()
print("Saved: figures/fig9_km_per_cancer.png")

# Save merged dataset
merged.to_csv("results/survival_dataset.tsv", sep="\t", index=False)
print(f"\nSaved: results/survival_dataset.tsv")
print("\nWeek 3 survival analysis complete!")
