import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# ── 1. Wood density dictionary: species -> (mean g/cm³, CV) ────────────────
DENSITY_DICT = {
    "Cassia fistula":       (0.96, 0.15),
    "Pongamia pinnata":     (0.64, 0.15),
    "Swietenia mahagoni":   (0.51, 0.15),
    "Terminalia bellirica": (0.76, 0.15),
}

# ── 2. Uncertainty constants ───────────────────────────────────────────────
DBH_CV         = 0.025   # 2.5 % normal
HEIGHT_CV      = 0.10    # 10 % normal (Changed from 0.15 as requested)
CHAVE_SIGMA_LN = 0.359   # Chave (2014) Eq.4 residual, lognormal multiplicative

# ── 3. Helpers ─────────────────────────────────────────────────────────────
def parse_dbh(val):
    if val is None:
        return []
    s = str(val).strip()
    if s == "" or s.lower() in ("nan", "na", "none"):
        return []
    matches = re.findall(r"=\s*(\d+\.?\d*)", s)
    return [float(n) for n in matches if float(n) > 0]

def chave_agb_e4(dbh_cm, rho_gcm3, height_m):
    D   = np.asarray(dbh_cm,  dtype=float)
    rho = np.asarray(rho_gcm3, dtype=float)
    H   = np.asarray(height_m, dtype=float)
    return 0.0673 * (rho * (D**2) * H)**0.976

def lognormal_from_mean_cv(mean, cv, size, rng):
    sigma = np.sqrt(np.log(1.0 + cv**2))
    mu    = np.log(mean) - 0.5 * sigma**2
    return rng.lognormal(mean=mu, sigma=sigma, size=size)

def monte_carlo_agb_tree(dbh_cm_list, rho_mean, rho_cv, height_m,
                         n_iter=5000, dbh_cv=DBH_CV, height_cv=HEIGHT_CV,
                         chave_sigma_ln=CHAVE_SIGMA_LN, seed=42):
    rng = np.random.default_rng(seed)
    dbh = np.asarray(dbh_cm_list, dtype=float)
    if dbh.size == 0:
        return np.zeros(n_iter)

    dbh_sd = np.maximum(dbh * dbh_cv, 1e-6)
    dbh_s  = rng.normal(loc=dbh, scale=dbh_sd, size=(n_iter, dbh.size))
    dbh_s  = np.clip(dbh_s, 1e-3, None)

    h_mean = float(height_m)
    h_s    = rng.normal(loc=h_mean, scale=max(h_mean * height_cv, 1e-6),
                        size=n_iter)
    h_s    = np.clip(h_s, 1e-2, None)

    rho_s  = lognormal_from_mean_cv(rho_mean, rho_cv, size=(n_iter, 1), rng=rng)
    agb_s  = chave_agb_e4(dbh_s, rho_s, h_s[:, None]).sum(axis=1)

    eps    = rng.lognormal(mean=0.0, sigma=chave_sigma_ln, size=n_iter)
    return agb_s * eps

# ── 4. Excel loader ────────────────────────────────────────────────────────
def load_excel(path):
    df = pd.read_excel(path, sheet_name=0)
    df = df.dropna(how="all")
    df = df[df["TreeID"].notna()].copy()

    def to_bool(v):
        if isinstance(v, (bool, np.bool_)): return bool(v)
        if isinstance(v, (int, float, np.integer, np.floating)):
            return not pd.isna(v) and v != 0
        return str(v).strip().upper() == "TRUE"

    df["TreeExists"] = df["TreeExists"].map(to_bool)
    df["TreeSpecies"] = df["TreeSpecies"].astype(str).str.strip()
    df = df[df["TreeExists"]].copy()

    df["DBH_eval"] = np.nan
    for idx, row in df.iterrows():
        dbh_raw = str(row.get("DBH ( in cm)", "")).strip()
        if dbh_raw.startswith("="):
            match = re.search(r"=([A-Z]+)(\d+)/PI\(\)", dbh_raw)
            if match:
                excel_row = int(match.group(2))
                pandas_idx = excel_row - 2
                if pandas_idx in df.index:
                    cbh_val = df.loc[pandas_idx, "CBH(cm)"]
                    if pd.notna(cbh_val):
                        df.at[idx, "DBH_eval"] = float(cbh_val) / np.pi
    return df

# ── 5. Pipeline (height now REQUIRED) ──────────────────────────────────────
def run_pipeline(df, n_iter=5000, seed=42):
    rows = []
    skipped = []

    for _, row in df.iterrows():
        tree_id = int(row["TreeID"])
        species = row["TreeSpecies"]

        if species not in DENSITY_DICT:
            skipped.append((tree_id, f"species not in DENSITY_DICT: {species}"))
            continue

        # ---- Height: now required, no fallback ----
        h_val = row.get("Height( in m)")
        if (pd.isna(h_val)
                or not isinstance(h_val, (int, float, np.integer, np.floating))
                or float(h_val) <= 0):
            skipped.append((tree_id, f"missing/invalid height ({h_val!r})"))
            continue
        height_m = float(h_val)

        # ---- DBH ----
        dbh_list = []
        stem_count = row.get("Stem Count", 1)
        if pd.notna(stem_count) and float(stem_count) > 1:
            dbh_list = parse_dbh(row.get("DBH ( in cm)"))
        else:
            if pd.notna(row["DBH_eval"]):
                dbh_list = [row["DBH_eval"]]
            else:
                cbh = row.get("CBH(cm)")
                if pd.notna(cbh) and isinstance(cbh, (int, float, np.integer, np.floating)):
                    dbh_list = [float(cbh) / np.pi]

        if not dbh_list:
            skipped.append((tree_id, "missing DBH/CBH"))
            continue

        rho_mean, rho_cv = DENSITY_DICT[species]
        samples = monte_carlo_agb_tree(
            dbh_list, rho_mean, rho_cv, height_m,
            n_iter=n_iter, seed=seed
        )

        rows.append({
            "TreeID":       tree_id,
            "Species":      species,
            "N_stems":      len(dbh_list),
            "DBH_mean_cm":  float(np.mean(dbh_list)),
            "Height_m":     height_m,
            "AGB_mean_kg":  float(samples.mean()),
            "AGB_sd_kg":    float(samples.std(ddof=1)),
            "AGB_p2.5_kg":  float(np.percentile(samples, 2.5)),
            "AGB_p97.5_kg": float(np.percentile(samples, 97.5)),
            "AGB_samples":  samples,
        })

    if skipped:
        print("\n=== Skipped trees ===")
        for tid, reason in skipped:
            print(f"[SKIP] TreeID {tid}: {reason}")

    return pd.DataFrame(rows)

# ── NEW: Export detailed stats to Excel ────────────────────────────────────
def export_to_excel(results_df, filename="agb_detailed_results.xlsx"):
    rows = []
    for _, row in results_df.iterrows():
        samples = row["AGB_samples"]
        
        # Calculate statistics
        recorded_agb = np.mean(samples) # The central deterministic / mean value
        median_agb   = np.median(samples)
        p2_5         = np.percentile(samples, 2.5)
        p97_5        = np.percentile(samples, 97.5)
        
        # CV = (Standard Deviation / Mean) * 100
        mean = samples.mean()
        sd   = samples.std(ddof=1)
        cv   = (sd / mean) * 100 if mean != 0 else 0.0

        rows.append({
            "TreeID":          row["TreeID"],
            "Species":         row["Species"],
            "Recorded_AGB_kg": round(recorded_agb, 3),
            "Mean_AGB_kg":     round(mean, 3),
            "Median_AGB_kg":   round(median_agb, 3),
            "Percentile_2.5":  round(p2_5, 3),
            "Percentile_97.5": round(p97_5, 3),
            "CV_%":            round(cv, 2)
        })
    
    df_out = pd.DataFrame(rows)
    df_out.to_excel(filename, index=False)
    print(f"Detailed stats saved to: {filename}")

# ── 6. Plots ───────────────────────────────────────────────────────────────
def make_plots(results_df, out_prefix="agb"):
    out_prefix = out_prefix.strip() 
    
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300,
        "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.spines.top": False, "axes.spines.right": False,
    })

    fig, ax = plt.subplots(figsize=(9, 4))
    x     = np.arange(len(results_df))
    means = results_df["AGB_mean_kg"].values
    lo    = means - results_df["AGB_p2.5_kg"].values
    hi    = results_df["AGB_p97.5_kg"].values - means
    ax.bar(x, means, yerr=[lo, hi], capsize=4,
           color="#4C72B0", edgecolor="black", linewidth=0.6,
           error_kw={"elinewidth": 1, "ecolor": "black"})
    ax.set_xticks(x)
    ax.set_xticklabels(results_df["TreeID"].astype(int).astype(str), rotation=90)
    ax.set_xlabel("Tree ID")
    ax.set_ylabel("Aboveground biomass (kg)")
    ax.set_title("Per-tree AGB (mean ± 95 % Monte Carlo CI)")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_per_tree_ci.png")
    plt.close(fig)

    n = len(results_df)
    n_cols = min(4, n)
    n_rows = int(np.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3 * n_cols, 2.6 * n_rows), squeeze=False)
    for i, (_, row) in enumerate(results_df.iterrows()):
        ax = axes[i // n_cols, i % n_cols]
        ax.hist(row["AGB_samples"], bins=40, color="#55A868", edgecolor="white", alpha=0.85)
        ax.axvline(row["AGB_mean_kg"], color="black", lw=1.2, ls="--",
                   label=f"mean = {row['AGB_mean_kg']:.2f} kg")
        ax.set_title(f"Tree {int(row['TreeID'])} — {row['Species']}", fontsize=8)
        ax.set_xlabel("AGB (kg)", fontsize=8)
        ax.set_ylabel("Count", fontsize=8)
        ax.legend(fontsize=7, loc="upper right")
    for j in range(n, n_rows * n_cols):
        axes[j // n_cols, j % n_cols].set_visible(False)
    fig.suptitle("Monte Carlo AGB distributions per tree", y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_distributions.png")
    plt.close(fig)

    groups, labels = [], []
    for _, row in results_df.iterrows():
        a = np.asarray(row["AGB_samples"], dtype=float).ravel()
        a = a[np.isfinite(a)]
        if a.size:
            groups.append(a)
            labels.append(str(int(row["TreeID"])))

    if not groups:
        print("[WARN] No valid samples for boxplot — skipping.")
        return

    fig, ax = plt.subplots(figsize=(max(8, 0.35 * len(groups) + 2), 4.5))
    bp = ax.boxplot(groups, patch_artist=True, widths=0.6, manage_ticks=False,
                    medianprops={"color": "black"},
                    showfliers=False, # ADDED: Hides the outliers from view
                    flierprops={"marker": "o", "markersize": 2, "alpha": 0.4})
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    colors = plt.colormaps["Set3"].resampled(len(groups))(np.arange(len(groups)))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_edgecolor("black")
    ax.set_xlabel("Tree ID")
    ax.set_ylabel("Aboveground biomass (kg)")
    ax.set_title("Per-tree AGB uncertainty — Monte Carlo (Chave et al., 2014)")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_per_tree_boxplot.png")
    plt.close(fig)

    print(f"Figures saved: {out_prefix}_per_tree_ci.png, "
          f"{out_prefix}_distributions.png, {out_prefix}_per_tree_boxplot.png")

# ── 7. Main ────────────────────────────────────────────────────────────────
def main():
    excel_path = "Eco_Report_Part1_v1 (1).xlsx"
    df = load_excel(excel_path)

    results = run_pipeline(df, n_iter=5000, seed=42)
    if results.empty:
        print("No trees with valid height + DBH — nothing to plot.")
        return

    # Save the basic summary CSV
    summary = results.drop(columns=["AGB_samples"]).round(3)
    print("\n=== AGB summary (Chave et al., 2014; Monte Carlo) ===")
    print(summary.to_string(index=False))
    summary.to_csv("agb_summary.csv", index=False)

    # NEW: Save the detailed stats Excel file
    export_to_excel(results, filename="agb_detailed_results.xlsx")

    # Generate plots
    make_plots(results, out_prefix="agb")

if __name__ == "__main__":
    main()
