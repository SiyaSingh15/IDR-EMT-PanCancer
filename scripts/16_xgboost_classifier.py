import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix
import xgboost as xgb
import shap

RESULTS_DIR = "results"
FIGURES_DIR = "figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

# Load master dataset
print("Loading master dataset...")
df = pd.read_csv("results/master_dataset.tsv", sep="\t")
print(f"Total samples: {len(df):,}")
print(f"EMT state distribution:\n{df['emt_state'].value_counts()}\n")

# Build feature matrix
# IDR mutation features - one binary flag per gene
idr_features = [
    "idr_mut_SNAI1","idr_mut_SNAI2","idr_mut_ZEB1","idr_mut_ZEB2",
    "idr_mut_TWIST1","idr_mut_TWIST2","idr_mut_ESRP1","idr_mut_ESRP2",
    "idr_mut_CDH1","idr_mut_VIM"
]

# Any EMT mutation flag
other_features = ["has_any_emt_mutation"]

# Cancer type one-hot encoding
ct_dummies = pd.get_dummies(df["cancer_type"], prefix="ct")

# Assemble feature matrix
X = pd.concat([
    df[idr_features].astype(int),
    df[other_features].astype(int),
    ct_dummies
], axis=1)

# Target: EMT state
y = df["emt_state"]
le = LabelEncoder()
y_enc = le.fit_transform(y)

print(f"Feature matrix: {X.shape[0]} samples x {X.shape[1]} features")
print(f"Class distribution: {dict(zip(le.classes_, np.bincount(y_enc)))}")
print(f"Feature names (first 15): {list(X.columns[:15])}")

# Train XGBoost with 5-fold stratified cross-validation
print("\nTraining XGBoost classifier (5-fold CV)...")

xgb_clf = xgb.XGBClassifier(
    n_estimators=200,
    max_depth=4,
    learning_rate=0.1,
    subsample=0.8,
    colsample_bytree=0.8,
    use_label_encoder=False,
    eval_metric="mlogloss",
    random_state=42,
    n_jobs=-1
)

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_scores = cross_val_score(xgb_clf, X, y_enc, cv=skf,
                            scoring="balanced_accuracy", n_jobs=-1)

print(f"\n5-fold CV balanced accuracy: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
print(f"Per-fold scores: {[f'{s:.3f}' for s in cv_scores]}")

# Fit on full dataset for SHAP
print("\nFitting on full dataset for SHAP analysis...")
xgb_clf.fit(X, y_enc)

# Full dataset predictions
y_pred = xgb_clf.predict(X)
print("\nClassification report (training set):")
print(classification_report(y_enc, y_pred, target_names=le.classes_))

# SHAP analysis
print("\nComputing SHAP values...")
explainer = shap.TreeExplainer(xgb_clf)
shap_values = explainer.shap_values(X)

# shap_values shape: (n_samples, n_features, n_classes) for multiclass
print(f"SHAP values shape: {np.array(shap_values).shape}")

# Mean absolute SHAP per feature per class
feature_names = list(X.columns)
n_classes = len(le.classes_)

# For each class, get mean |SHAP|
shap_importance = pd.DataFrame()
for i, cls in enumerate(le.classes_):
    if isinstance(shap_values, list):
        sv = shap_values[i]
    else:
        sv = shap_values[:, :, i]
    mean_abs = np.abs(sv).mean(axis=0)
    shap_importance[cls] = mean_abs

shap_importance.index = feature_names
shap_importance["mean_all"] = shap_importance.mean(axis=1)
shap_importance = shap_importance.sort_values("mean_all", ascending=False)

print("\nTop 15 features by mean |SHAP|:")
print(shap_importance[["Epithelial","Hybrid","Mesenchymal","mean_all"]].head(15).to_string())

# Figure E: SHAP feature importance - IDR features only
print("\nGenerating Figure E: SHAP importance for IDR features...")

# Focus on IDR mutation features
idr_shap = shap_importance.loc[
    [f for f in idr_features if f in shap_importance.index]
].copy()

# Clean feature names for plot
idr_shap.index = [f.replace("idr_mut_","") for f in idr_shap.index]
idr_shap = idr_shap.sort_values("mean_all", ascending=True)

colors_cls = {
    "Epithelial": "#3498DB",
    "Hybrid": "#F39C12",
    "Mesenchymal": "#E74C3C"
}

fig, ax = plt.subplots(figsize=(9, 6))
y_pos = np.arange(len(idr_shap))
width = 0.25

for i, cls in enumerate(le.classes_):
    ax.barh(y_pos + i*width - width, idr_shap[cls],
            width, label=cls, color=colors_cls[cls], alpha=0.85)

ax.set_yticks(y_pos)
ax.set_yticklabels(idr_shap.index, fontsize=11)
ax.set_xlabel("Mean |SHAP value|", fontsize=11)
ax.set_title("SHAP feature importance: IDR mutation features\n"
             "by predicted EMT state (XGBoost, 5-fold CV)", fontsize=11)
ax.legend(title="Predicted class", fontsize=9)
ax.axvline(0, color="black", linewidth=0.5)
plt.tight_layout()
plt.savefig("figures/figE_shap_idr_importance.png", dpi=150)
plt.close()
print("Saved: figures/figE_shap_idr_importance.png")

# Figure F: SHAP beeswarm-style plot - top 20 features overall
print("Generating Figure F: SHAP summary plot...")

top20 = shap_importance.head(20).index.tolist()
X_top20 = X[top20]

# Use class index 2 (Mesenchymal) for beeswarm 
if isinstance(shap_values, list):
    sv_mesen = shap_values[list(le.classes_).index("Mesenchymal")]
else:
    sv_mesen = shap_values[:, :, list(le.classes_).index("Mesenchymal")]

sv_top20 = sv_mesen[:, [list(X.columns).index(f) for f in top20]]

fig, ax = plt.subplots(figsize=(9, 8))
shap.summary_plot(
    sv_top20, X_top20,
    feature_names=[f.replace("idr_mut_","IDR:").replace("ct_","").replace("has_any_emt_mutation","any EMT mut")
                   for f in top20],
    plot_type="dot",
    show=False,
    max_display=20
)
plt.title("SHAP beeswarm: top 20 features predicting Mesenchymal state", fontsize=11)
plt.tight_layout()
plt.savefig("figures/figF_shap_beeswarm.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved: figures/figF_shap_beeswarm.png")

# Figure G: Confusion matrix
print("Generating Figure G: confusion matrix...")

cm = confusion_matrix(y_enc, y_pred)
fig, ax = plt.subplots(figsize=(6, 5))
im = ax.imshow(cm, cmap="Blues")
plt.colorbar(im, ax=ax)
ax.set_xticks(range(n_classes))
ax.set_yticks(range(n_classes))
ax.set_xticklabels(le.classes_, fontsize=10)
ax.set_yticklabels(le.classes_, fontsize=10)
ax.set_xlabel("Predicted", fontsize=11)
ax.set_ylabel("True", fontsize=11)
ax.set_title(f"XGBoost confusion matrix\n(5-fold CV balanced accuracy = {cv_scores.mean():.3f})", fontsize=11)
for i in range(n_classes):
    for j in range(n_classes):
        ax.text(j, i, str(cm[i,j]), ha="center", va="center",
                color="white" if cm[i,j] > cm.max()/2 else "black", fontsize=12)
plt.tight_layout()
plt.savefig("figures/figG_confusion_matrix.png", dpi=150)
plt.close()
print("Saved: figures/figG_confusion_matrix.png")

# Figure H: CV accuracy per fold
fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(range(1, 6), cv_scores, color="#2980B9", alpha=0.8, width=0.6)
ax.axhline(cv_scores.mean(), color="#E74C3C", linestyle="--",
           linewidth=2, label=f"Mean = {cv_scores.mean():.3f}")
ax.axhline(1/3, color="#AAAAAA", linestyle=":", linewidth=1.5,
           label="Random baseline = 0.333")
ax.set_xlabel("CV fold", fontsize=11)
ax.set_ylabel("Balanced accuracy", fontsize=11)
ax.set_title("XGBoost 5-fold cross-validation performance\n"
             "IDR mutation + cancer type → E/Hybrid/M prediction", fontsize=11)
ax.set_ylim(0, 1)
ax.legend(fontsize=9)
ax.set_xticks(range(1, 6))
plt.tight_layout()
plt.savefig("figures/figH_cv_performance.png", dpi=150)
plt.close()
print("Saved: figures/figH_cv_performance.png")

# Save model results
results = {
    "cv_balanced_accuracy_mean": float(cv_scores.mean()),
    "cv_balanced_accuracy_std": float(cv_scores.std()),
    "cv_scores_per_fold": cv_scores.tolist(),
    "n_features": int(X.shape[1]),
    "n_samples": int(X.shape[0]),
    "top_idr_features_by_shap": idr_shap.sort_values("mean_all",ascending=False).head(5)["mean_all"].to_dict()
}

import json
with open("results/xgboost_results.json", "w") as f:
    json.dump(results, f, indent=2)

shap_importance.to_csv("results/shap_importance.tsv", sep="\t")

print(f"\n=== XGBOOST SUMMARY ===")
print(f"5-fold CV balanced accuracy: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
print(f"Random baseline: 0.333")
print(f"Improvement over random: {(cv_scores.mean()-0.333)/0.333*100:.1f}%")
print(f"\nTop IDR features by SHAP (Mesenchymal prediction):")
print(idr_shap.sort_values("mean_all",ascending=False)[["Mesenchymal","mean_all"]].head(5).to_string())
print("\nAll figures saved to figures/")
print("Results saved to results/xgboost_results.json")
