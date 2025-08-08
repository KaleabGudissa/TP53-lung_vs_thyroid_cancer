# TP53 in LUAD vs THCA (TCGA PanCancer Atlas) — Single‑Folder Repo

Compare **TP53 mutation frequency** in Lung Adenocarcinoma (**LUAD**) vs Thyroid Carcinoma (**THCA**), and test whether **TP53 status associates with overall survival in LUAD**. The repo uses a **single‑folder layout** — all inputs and outputs live in the repo root.

---

# 🚀 Quick Start (flat repo)

## 1) Put the study bundles in the repo root

Download and Place these files in data/raw :

* `luad_tcga_pan_can_atlas_2018.tar.gz` from https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018 
* `thca_tcga_pan_can_atlas_2018.tar.gz` from https://www.cbioportal.org/study/summary?id=thca_tcga_pan_can_atlas_2018


## 2) Create an environment and install requirements

```bash
python -m venv .venv
# mac/linux
source .venv/bin/activate
# windows (powershell)
# .venv\Scripts\Activate.ps1

python -m pip install --upgrade pip
pip install -r requirements.txt
```

## 3) Run the pipeline

```bash
python run_all.py
```



# 📄 Outputs (created in repo root)

## Figures

* `fig1_tp53_pct_luad_vs_thca.png` and `fig1_tp53_pct_luad_vs_thca.pdf` — % of patients with non‑silent **TP53** mutations in LUAD vs THCA with 95% **Wilson CIs**.
* `fig2_km_luad_tp53.png` and `fig2_km_luad_tp53.pdf` — **Kaplan–Meier** overall survival in LUAD, TP53‑mutated vs wild‑type, with **log‑rank p‑value**.

## Tables (CSV)

* `table1_tp53_counts.csv` — cohort size, TP53+, percentage, Wilson 95% CI.
* `table2_tp53_stats.csv` — Chi‑square and Fisher’s exact tests, **Odds Ratio (OR)**, **Risk Ratio (RR)** with 95% CIs.
* `table3_cox_summary.csv` — (best‑effort) **Cox PH** summary (unadjusted TP53 HR; adds age/sex/stage if present).

## Supplementary / processed (for transparency)

* `luad_tp53_status.csv`, `thca_tp53_status.csv` — per‑patient TP53 status (1 if any non‑silent TP53 variant, else 0).
* `luad_tp53_rows.csv`, `thca_tp53_rows.csv` — all TP53 mutation rows from each study.
* `luad_data_mutations.txt`, `thca_data_mutations.txt` — extracted mutation tables analyzed.
* `luad_data_clinical_patient.txt`, `thca_data_clinical_patient.txt` — extracted patient clinical tables.

---

# 🧪 Methods (summary)

* **Source:** cBioPortal TCGA PanCancer Atlas study bundles (LUAD & THCA).
* **Patient ID:** first 12 chars of `Tumor_Sample_Barcode`.
* **TP53 mutated:** any **non‑silent** `Variant_Classification` (Missense, Nonsense, Frame\_Shift\_Ins/Del, Splice\_Site, Translation\_Start\_Site, Nonstop, In\_Frame\_Ins/Del).
* **Frequency:** % TP53‑mutated with **Wilson 95% CIs**.
* **Between‑cancer test:** Chi‑square (no Yates) and Fisher’s exact (for small THCA counts). Report **OR** and **RR** with 95% CIs.
* **Survival (LUAD):** Kaplan–Meier + log‑rank; optional **Cox PH** (unadjusted; includes age/sex/stage if present).

---

# ✅ Expected headline numbers (sanity check)

* **Mutation‑frequency cohort sizes:** LUAD ≈ 562, THCA ≈ 484 (total ≈ 1,046).
* **TP53 mutated:** LUAD ≈ **51%** (≈287/562), THCA ≈ **0.4%** (≈2/484).
* **LUAD survival:** log‑rank **p ≈ 0.13** (typically not significant). *(Your exact CSV outputs are the source of truth.)*

---

# 🧹 Don’t commit raw bundles

GitHub blocks files **> 100 MB**. Keep large files out of git.

## Minimal `.gitignore`

```gitignore
*.tar.gz
*.maf
*.gz
```

If you accidentally committed them, remove from history (e.g., `git filter-repo`) and force‑push; keep the raw files local only.

---

# 🔧 Troubleshooting

* **Push blocked (>100MB):** you committed a `.tar.gz`. Untrack it and rewrite history; keep `.tar.gz` out of git.
* **Module not found:** run `pip install -r requirements.txt`.
* **KM/Cox missing:** if one LUAD group has zero OS entries, the script will skip KM/Cox but still write frequency tables/figures.

---

# ℹ️ Citation 

* **cBioPortal for Cancer Genomics** — TCGA PanCancer Atlas studies.  https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018 
&  https://www.cbioportal.org/study/summary?id=thca_tcga_pan_can_atlas_2018
* **pandas, numpy, matplotlib, SciPy, lifelines** — analysis, statistics, and plots. 
