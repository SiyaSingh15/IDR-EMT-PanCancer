# IDR Mutations in EMT Master Regulators as Modulators of Phenotypic Plasticity

A pan-cancer TCGA computational study

[![Python](https://img.shields.io/badge/Python-3.10-blue)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

First systematic pan-cancer analysis of somatic mutations within intrinsically
disordered regions (IDRs) of 10 EMT regulator proteins across 10,295 TCGA
tumour samples.

**Proteins studied:** SNAI1, SNAI2, ZEB1, ZEB2, TWIST1, TWIST2, ESRP1, ESRP2, CDH1, VIM

## Key Findings

- 255 IDR-resident somatic mutations across 222 unique patients (15.3% of coding EMT mutations)
- ESRP1 IDR mutations deplete Hybrid E/M state (28% vs 50% WT; Fisher's p=0.096)
- TCGA-UCEC IDR-mutant patients show significantly better survival (log-rank p=0.007)
- Mesenchymal EMT state predicts worse survival independently (Cox HR=1.20, p<0.005)
- XGBoost + SHAP: ZEB1/ZEB2 IDR mutations predict Hybrid state; ESRP1/ESRP2 predict Mesenchymal

## Repository Structure
```
IDR-EMT-PanCancer/
├── scripts/          # Analysis pipeline (scripts 05-16)
├── results/          # Key output tables and statistics
├── figures/          # All publication figures
├── requirements.txt  # Python dependencies
└── README.md
```

## Pipeline

| Script | Function |
|--------|----------|
| `05_build_mutation_table.py` | Load MC3 MAF, filter EMT genes |
| `06_fetch_idr_boundaries.py` | Query MobiDB API for IDR coordinates |
| `07_flag_idr_mutations.py` | Map mutations onto IDR/ordered domain |
| `08_fetch_expression_cbioportal.py` | Fetch RNA-seq from cBioPortal API |
| `09_compute_emt_scores.py` | Compute EMP scores, classify E/Hybrid/M |
| `10_build_master_dataset.py` | Merge mutations + EMT scores |
| `11_statistical_analysis.py` | Pan-cancer figures + Mann-Whitney |
| `12_stratified_analysis.py` | Within-cancer + per-gene tests + FDR |
| `13_fix_cancer_types.py` | GDC barcode-to-cancer-type mapping |
| `14_survival_analysis.py` | KM + log-rank + Cox PH |
| `15_final_figures.py` | Publication figure generation |
| `16_xgboost_classifier.py` | XGBoost + SHAP analysis |

## Data Availability

Data files are not tracked due to size.

**Mutation data (MC3):**
```
https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc
```
Save to `data/mafs/mc3.v0.2.8.PUBLIC.maf.gz`

**Expression data:** Auto-retrieved via cBioPortal API (`script 08`)

**IDR annotations:** Auto-retrieved via MobiDB REST API (`script 06`)

## Installation
```bash
git clone https://github.com/SiyaSingh15/IDR-EMT-PanCancer.git
cd IDR-EMT-PanCancer
conda create -n idr_emt python=3.10
conda activate idr_emt
pip install -r requirements.txt
```


## Author

**Siya** | IISER Tirupati | April 2026

## License

MIT — see [LICENSE](LICENSE)
