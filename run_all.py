#!/usr/bin/env python3
import tarfile, re, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from scipy.stats import chi2_contingency, fisher_exact

RAW = Path("data/raw")
PROC = Path("data/processed")
FIGS = Path("results/figures")
TABLES = Path("results/tables")
DOCS = Path("docs")
for d in [PROC, FIGS, TABLES, DOCS]:
    d.mkdir(parents=True, exist_ok=True)

NON_SILENT = {
    "Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins",
    "Splice_Site","Translation_Start_Site","Nonstop_Mutation","In_Frame_Del","In_Frame_Ins"
}

def find_tar(name_hint: str) -> Path:
    hits = sorted([p for p in RAW.glob("*.tar.gz") if name_hint.lower() in p.name.lower()])
    if not hits:
        raise FileNotFoundError(f"Put a *{name_hint}*.tar.gz in {RAW} (from cBioPortal Datasets).")
    return hits[0]

def extract_member(tar_path: Path, patterns, out_path: Path) -> Path:
    with tarfile.open(tar_path, "r:gz") as tar:
        member = None
        for patt in patterns:
            for m in tar.getmembers():
                if re.search(patt, m.name):
                    member = m; break
            if member: break
        if member is None:
            raise FileNotFoundError(f"No member matching {patterns} in {tar_path.name}")
        with tar.extractfile(member) as fsrc:
            out_path.write_bytes(fsrc.read())
    return out_path

def load_mutations(path: Path) -> pd.DataFrame:
    with open(path, "rt", errors="ignore") as f:
        header = f.readline().rstrip("\n").split("\t")
    need = ["Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification"]
    use = [c for c in need if c in header]
    if len(use) < 3:
        raise ValueError(f"{path} missing required columns {need}")
    df = pd.read_csv(path, sep="\t", usecols=use, dtype=str, low_memory=False)
    df["patient"] = df["Tumor_Sample_Barcode"].astype(str).str[:12]
    return df

def load_cbio_patient_clin(path: Path) -> pd.DataFrame:
    rows = []
    with open(path, "rt", errors="ignore") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            rows.append(line.rstrip("\n").split("\t"))
    if not rows:
        raise ValueError(f"{path} appears empty.")
    header = rows[0]; data = rows[1:]
    df = pd.DataFrame(data, columns=header)

    def to_event(x):
        if pd.isna(x): return pd.NA
        v = str(x).upper()
        if "DECEASED" in v or v.startswith("1"): return 1
        if "LIVING" in v or v.startswith("0"): return 0
        return pd.NA

    out = pd.DataFrame({
        "patient_id": df.get("PATIENT_ID"),
        "os_months": pd.to_numeric(df.get("OS_MONTHS"), errors="coerce"),
        "os_event": df.get("OS_STATUS").apply(to_event)
    })
    # optional covariates if present
    if "AGE" in df.columns:
        out["age"] = pd.to_numeric(df["AGE"], errors="coerce")
    if "SEX" in df.columns:
        out["sex"] = df["SEX"]
    if "AJCC_PATHOLOGIC_STAGE" in df.columns:
        out["stage"] = df["AJCC_PATHOLOGIC_STAGE"]
    elif "PATHOLOGIC_STAGE" in df.columns:
        out["stage"] = df["PATHOLOGIC_STAGE"]
    return out

def tp53_status(mut_df: pd.DataFrame):
    tp = mut_df[(mut_df["Hugo_Symbol"]=="TP53") & (mut_df["Variant_Classification"].isin(NON_SILENT))].copy()
    mutated = set(tp["patient"])
    pats = mut_df["patient"].dropna().unique()
    status = pd.DataFrame({"patient_id": pats})
    status["tp53_mut"] = status["patient_id"].apply(lambda x: 1 if x in mutated else 0)
    return status, tp

def wilson_ci(k, n):
    if n == 0: return (np.nan, np.nan)
    z = 1.959963984540054
    p = k / n
    denom = 1 + z**2/n
    center = (p + z*z/(2*n)) / denom
    half = (z * math.sqrt((p*(1-p)/n) + (z*z/(4*n*n)))) / denom
    return (center - half, center + half)

def odds_ratio_ci(a,b,c,d):
    if min(a,b,c,d) == 0:
        a+=0.5; b+=0.5; c+=0.5; d+=0.5
    or_val = (a*d)/(b*c)
    se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    z = 1.959963984540054
    lo = math.exp(math.log(or_val) - z*se)
    hi = math.exp(math.log(or_val) + z*se)
    return or_val, lo, hi

def risk_ratio_ci(a,b,c,d):
    if (a+b)==0 or (c+d)==0:
        return (np.nan, np.nan, np.nan)
    p1 = a/(a+b); p0 = c/(c+d)
    rr = p1/p0 if p0>0 else np.inf
    if min(a,b,c,d) == 0:
        a+=0.5; b+=0.5; c+=0.5; d+=0.5
    se = math.sqrt( (1/a) - (1/(a+b)) + (1/c) - (1/(c+d)) )
    z = 1.959963984540054
    lo = math.exp(math.log(rr) - z*se)
    hi = math.exp(math.log(rr) + z*se)
    return rr, lo, hi

def main():
    # 1) Extract needed files
    luad_tar = find_tar("luad")
    thca_tar = find_tar("thca")

    luad_mut_p = extract_member(luad_tar, [r"/data_mutations_extended\.txt$", r"/data_mutations\.txt$"], PROC/"luad_data_mutations.txt")
    thca_mut_p = extract_member(thca_tar, [r"/data_mutations_extended\.txt$", r"/data_mutations\.txt$"], PROC/"thca_data_mutations.txt")

    luad_clin_p = extract_member(luad_tar, [r"/data_clinical_patient\.txt$", r"/data_clinical\.txt$"], PROC/"luad_data_clinical_patient.txt")
    thca_clin_p = extract_member(thca_tar, [r"/data_clinical_patient\.txt$", r"/data_clinical\.txt$"], PROC/"thca_data_clinical_patient.txt")

    # 2) TP53 status per patient
    luad_mut = load_mutations(luad_mut_p)
    thca_mut = load_mutations(thca_mut_p)
    luad_status, luad_tp53_rows = tp53_status(luad_mut)
    thca_status, thca_tp53_rows = tp53_status(thca_mut)

    # Save supplements
    luad_status.to_csv(TABLES/"luad_tp53_status.csv", index=False)
    thca_status.to_csv(TABLES/"thca_tp53_status.csv", index=False)
    luad_tp53_rows.to_csv(TABLES/"luad_tp53_rows.csv", index=False)
    thca_tp53_rows.to_csv(TABLES/"thca_tp53_rows.csv", index=False)

    # 3) Table 1 — counts/percents + Wilson CI
    luad_n = luad_status["patient_id"].nunique()
    luad_m = int(luad_status["tp53_mut"].sum())
    thca_n = thca_status["patient_id"].nunique()
    thca_m = int(thca_status["tp53_mut"].sum())

    luad_ci = wilson_ci(luad_m, luad_n)
    thca_ci = wilson_ci(thca_m, thca_n)

    table1 = pd.DataFrame({
        "Cohort": ["LUAD","THCA"],
        "Patients with mutation data (n)": [luad_n, thca_n],
        "Patients with non-silent TP53 (n)": [luad_m, thca_m],
        "Percent TP53 mutated (%)": [round(100*luad_m/luad_n,1), round(100*thca_m/thca_n,1)],
        "Wilson 95% CI (lower %)": [round(100*luad_ci[0],1), round(100*thca_ci[0],3)],
        "Wilson 95% CI (upper %)": [round(100*luad_ci[1],1), round(100*thca_ci[1],3)],
    })
    table1.to_csv(TABLES/"table1_tp53_counts.csv", index=False)

    # 4) Frequency tests — chi-square + Fisher + effect sizes
    a = luad_m; b = luad_n - luad_m
    c = thca_m; d = thca_n - thca_m
    table = np.array([[a,b],[c,d]])

    chi2, chi2_p, _, _ = chi2_contingency(table, correction=False)
    fisher_or, fisher_p = (np.nan, np.nan)
    if (table < 5).any():
        fisher_or, fisher_p = fisher_exact(table, alternative="two-sided")

    or_val, or_lo, or_hi = odds_ratio_ci(a,b,c,d)
    rr_val, rr_lo, rr_hi = risk_ratio_ci(a,b,c,d)

    table2 = pd.DataFrame({
        "Stat": ["Chi-square (no Yates)", "Chi-square p-value",
                 "Fisher’s exact OR", "Fisher’s exact p-value",
                 "Odds Ratio (Wald CI)", "OR 95% CI lower", "OR 95% CI upper",
                 "Risk Ratio (Katz CI)", "RR 95% CI lower", "RR 95% CI upper"],
        "Value": [round(chi2,4), chi2_p,
                  fisher_or, fisher_p,
                  or_val, or_lo, or_hi,
                  rr_val, rr_lo, rr_hi]
    })
    table2.to_csv(TABLES/"table2_tp53_stats.csv", index=False)

    # 5) Figure 1 — % TP53 with 95% CI
    xs = ["LUAD","THCA"]
    props = [100*luad_m/luad_n, 100*thca_m/thca_n]
    lows = [100*(luad_m/luad_n - luad_ci[0]), 100*(thca_m/thca_n - thca_ci[0])]
    highs = [100*(luad_ci[1] - luad_m/luad_n), 100*(thca_ci[1] - thca_m/thca_n)]
    asym_err = np.array([lows, highs])

    plt.figure()
    plt.bar(xs, props, yerr=asym_err, capsize=6)
    for i, v in enumerate(props):
        plt.text(i, v, f"{v:.1f}%", ha="center", va="bottom")
    plt.ylabel("% patients with non-silent TP53")
    plt.title("TP53 mutation frequency: LUAD vs THCA (95% CI)")
    plt.tight_layout()
    plt.savefig(FIGS/"fig1_tp53_pct_luad_vs_thca.png", dpi=300)
    plt.savefig(FIGS/"fig1_tp53_pct_luad_vs_thca.pdf")
    plt.close()

    # 6) LUAD survival — KM + log-rank + Cox (with encoding)
    clin_luad = load_cbio_patient_clin(luad_clin_p)
    km = clin_luad.merge(luad_status, on="patient_id", how="inner")

    # Ensure numeric types for lifelines
    km["os_months"] = pd.to_numeric(km["os_months"], errors="coerce")
    km["os_event"]  = pd.to_numeric(km["os_event"],  errors="coerce")
    km = km.dropna(subset=["os_months","os_event"])
    km["os_months"] = km["os_months"].astype(float)
    km["os_event"]  = km["os_event"].astype(int)

    g1 = km[km["tp53_mut"]==1]
    g0 = km[km["tp53_mut"]==0]

    if len(g1) == 0 or len(g0) == 0:
        (TABLES/"table3_cox_summary.csv").write_text("KM/Cox skipped: one group empty.\n", encoding="utf-8")
    else:
        # KM + log-rank
        kmf0, kmf1 = KaplanMeierFitter(), KaplanMeierFitter()
        kmf0.fit(durations=g0["os_months"], event_observed=g0["os_event"], label=f"TP53 wild-type (n={len(g0)})")
        kmf1.fit(durations=g1["os_months"], event_observed=g1["os_event"], label=f"TP53 mutated (n={len(g1)})")
        res = logrank_test(g1["os_months"], g0["os_months"], event_observed_A=g1["os_event"], event_observed_B=g0["os_event"])
        chi2_lr, p_lr = float(res.test_statistic), float(res.p_value)

        ax = kmf0.plot(ci_show=False)
        kmf1.plot(ax=ax, ci_show=False)
        ax.set_xlabel("Overall survival (months)")
        ax.set_ylabel("Survival probability")
        ax.set_title(f"LUAD KM by TP53 (log-rank chi2={chi2_lr:.2f}, p={p_lr:.3g})")
        plt.tight_layout()
        plt.savefig(FIGS/"fig2_km_luad_tp53.png", dpi=300)
        plt.savefig(FIGS/"fig2_km_luad_tp53.pdf")
        plt.close()

        # Cox PH (encode categorical covariates to numeric and fit)
        try:
            cox_df = km[["os_months","os_event","tp53_mut"]].copy()

            # age if present
            if "age" in km.columns:
                cox_df["age"] = pd.to_numeric(km["age"], errors="coerce")

            # sex -> numeric (Female=0, Male=1)
            if "sex" in km.columns:
                cox_df["sex_num"] = km["sex"].astype(str).str.upper().map({"FEMALE":0, "MALE":1})

            # stage -> ordinal 1..4 (handles IA/IB/IIIA etc.)
            if "stage" in km.columns:
                def stage_to_num(s):
                    s = str(s).upper()
                    if "IV" in s:   return 4
                    if "III" in s:  return 3
                    if "II" in s:   return 2
                    if "I" in s:    return 1
                    return np.nan
                cox_df["stage_num"] = km["stage"].apply(stage_to_num)

            # keep only numeric columns and drop rows with NA
            num_cols = [c for c in cox_df.columns if c not in {"sex","stage"}]
            cox_df = cox_df[num_cols].apply(pd.to_numeric, errors="coerce").dropna()

            if len(cox_df) < 30:
                (TABLES/"table3_cox_summary.csv").write_text(f"Cox skipped: not enough complete rows ({len(cox_df)}).\n", encoding="utf-8")
            else:
                cph = CoxPHFitter()
                cph.fit(cox_df, duration_col="os_months", event_col="os_event")
                cph.summary.reset_index().to_csv(TABLES/"table3_cox_summary.csv", index=False)
        except Exception as e:
            (TABLES/"table3_cox_summary.csv").write_text(f"Could not fit Cox model: {e}\n", encoding="utf-8")

    # 7) Minimal provenance README (UTF-8)
    (DOCS/"DATASET_README.txt").write_text(
"""Inputs: cBioPortal TCGA PanCancer Atlas study bundles
- luad_tcga_pan_can_atlas_2018.tar.gz
- thca_tcga_pan_can_atlas_2018.tar.gz

Extracted tables analyzed:
- data_mutations.txt (or data_mutations_extended.txt)
- data_clinical_patient.txt

Definitions:
- patient_id = first 12 chars of Tumor_Sample_Barcode
- TP53 mutated = any non-silent Variant_Classification

Outputs:
- Table1: counts, %, Wilson 95% CI
- Table2: chi-square + Fisher, OR, RR (95% CI)
- Figure1: % TP53 (95% CI bars)
- Figure2: KM (log-rank)
- Table3: Cox PH (if possible)
""",
        encoding="utf-8"
    )

    # Console summary
    print("\n=== SUMMARY ===")
    print(f"LUAD: {luad_m}/{luad_n} TP53-mutated ({100*luad_m/luad_n:.1f}%)")
    print(f"THCA: {thca_m}/{thca_n} TP53-mutated ({100*thca_m/thca_n:.1f}%)")
    print(f"Chi-square: chi2={chi2:.2f}, p={chi2_p:.3g} | Fisher p={fisher_p if not np.isnan(fisher_p) else 'n/a'}")
    print(f"OR={or_val:.2f} [{or_lo:.2f}, {or_hi:.2f}] | RR={rr_val:.2f} [{rr_lo:.2f}, {rr_hi:.2f}]")
    if (FIGS/"fig2_km_luad_tp53.png").exists():
        print("KM/Cox completed (see figures/tables).")
    else:
        print("KM skipped (one group empty).")
    print("Outputs → results/figures, results/tables, docs/")
    print("=============\n")

if __name__ == "__main__":
    main()
