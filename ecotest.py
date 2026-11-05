import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ── 1. Wood density dictionary (g/cm³) ─────────────────────────────────────
DENSITY_DICT = {
    "Cassia fistula":        [0.71, 0.96],
    "Pongamia pinnata":      [0.54, 0.67],
    "Swietenia mahagoni":    [0.54, 0.83],
    "Terminalia bellirica":  [0.57, 0.90],
}

# ── 2. Helpers ─────────────────────────────────────────────────────────────
def parse_dbh(val):
    """Return list of positive floats found AFTER an equals sign."""
    if val is None:
        return []
    s = str(val).strip()
    if s == "" or s.lower() in ("nan", "na", "none"):
        return []
    
    # This regex looks for numbers immediately following an '=' sign.
    # It ignores the numbers in "DBH1", "DBH2", etc.
    matches = re.findall(r"=\s*(\d+\.?\d*)", s)
    
    return [float(n) for n in matches if float(n) > 0]

def chave_agb_e4(dbh_cm, rho_gcm3, height_m):
    """Chave et al. (2014) Eq. 4: 0.0673 * (ρ * D² * H)^0.976"""
    D   = np.asarray(dbh_cm,  dtype=float)
    rho = np.asarray(rho_gcm3, dtype=float)
    H   = np.asarray(height_m, dtype=float)
    return 0.0673 * (rho * (D**2) * H)**0.976

def monte_carlo_agb_tree(dbh_cm_list, rho_min, rho_max, height_m,
                         n_iter=5000, dbh_frac_err=0.02, seed=42):
    """Per-tree MC: DBH ±2%, density uniform, height fixed."""
    rng = np.random.default_rng(seed)
    dbh = np.asarray(dbh_cm_list, dtype=float)
    if dbh.size == 0:
        return np.zeros(n_iter)
    
    dbh_s = dbh * rng.uniform(1 - dbh_frac_err, 1 + dbh_frac_err,
                              size=(n_iter, dbh.size))
    rho_s = rng.uniform(rho_min, rho_max, size=(n_iter, dbh.size))
    
    return chave_agb_e4(dbh_s, rho_s, height_m).sum(axis=1)

# ── 3. Excel loader with Formula Evaluation ────────────────────────────────
def load_excel(path):
    df = pd.read_excel(path, sheet_name=0)
    df = df.dropna(how="all")
    df = df[df["TreeID"].notna()].copy()

    # Parse TreeExists (handles 1.0/0.0 or TRUE/FALSE)
    def to_bool(v):
        if isinstance(v, (bool, np.bool_)): return bool(v)
        if isinstance(v, (int, float, np.integer, np.floating)):
            return not pd.isna(v) and v != 0
        return str(v).strip().upper() == "TRUE"
        
    df["TreeExists"] = df["TreeExists"].map(to_bool)
    df["TreeSpecies"] = df["TreeSpecies"].astype(str).str.strip()
    df = df[df["TreeExists"]].copy()

    # Evaluate Excel formulas in the DBH column safely
    df["DBH_eval"] = np.nan
    for idx, row in df.iterrows():
        dbh_raw = str(row.get("DBH ( in cm)", "")).strip()
        if dbh_raw.startswith("="):
            # Excel formula like "=F2/PI()"
            match = re.search(r"=([A-Z]+)(\d+)/PI\(\)", dbh_raw)
            if match:
                excel_row = int(match.group(2))
                # Excel row 2 is pandas index 0, so subtract 2
                pandas_idx = excel_row - 2
                if pandas_idx in df.index:
                    cbh_val = df.loc[pandas_idx, "CBH(cm)"]
                    if pd.notna(cbh_val):
                        df.at[idx, "DBH_eval"] = float(cbh_val) / np.pi
                        
    return df

# ── 4. Pipeline ────────────────────────────────────────────────────────────
def run_pipeline(df, n_iter=5000, seed=42):
    rows = []
    for _, row in df.iterrows():
        species = row["TreeSpecies"]
        if species not in DENSITY_DICT:
            continue
            
        dbh_list = []
        stem_count = row.get("Stem Count", 1)
        
        if pd.notna(stem_count) and float(stem_count) > 1:
            # MULTI-STEM: parse strings like "DBH1=0.767 DBH2=0.318"
            dbh_list = parse_dbh(row.get("DBH ( in cm)"))
        else:
            # SINGLE-STEAM: use the evaluated formula or fall back to CBH/pi
            if pd.notna(row["DBH_eval"]):
                dbh_list = [row["DBH_eval"]]
            else:
                cbh = row.get("CBH(cm)")
                if pd.notna(cbh) and isinstance(cbh, (int, float)):
                    dbh_list = [float(cbh) / np.pi]
                    
        if not dbh_list:
            continue

        rho_min, rho_max = DENSITY_DICT[species]
        
        # Get Height (fallback to 5.0 meters if blank or multi-stem 1.66m)
        h_val = row.get("Height( in m)")
        if pd.isna(h_val) or not isinstance(h_val, (int, float)):
            height_m = 5.0
        else:
            height_m = float(h_val)

        samples = monte_carlo_agb_tree(dbh_list, rho_min, rho_max, height_m,
                                       n_iter=n_iter, seed=seed)
                                       
        rows.append({
            "TreeID":       int(row["TreeID"]),
            "Species":      species,
            "N_stems":      len(dbh_list),
            "DBH_mean_cm":  float(np.mean(dbh_list)),
            "AGB_mean_kg":  float(samples.mean()),
            "AGB_sd_kg":    float(samples.std(ddof=1)),
            "AGB_p2.5_kg":  float(np.percentile(samples, 2.5)),
            "AGB_p97.5_kg": float(np.percentile(samples, 97.5)),
            "AGB_samples":  samples,
        })
    return pd.DataFrame(rows)

# ── 5. Plots ───────────────────────────────────────────────────────────────
def make_plots(results_df, out_prefix="agb"):
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300,
        "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.spines.top": False, "axes.spines.right": False,
    })

    # Fig 1: per-tree mean AGB + 95 % MC CI
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

    # Fig 2: per-tree posterior histograms
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

    # Fig 3: PER-TREE boxplot of Monte Carlo AGB
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

    print(f"Figures saved: {out_prefix}_per_tree_ci.png, {out_prefix}_distributions.png, {out_prefix}_per_tree_boxplot.png")

# ── 6. Main ────────────────────────────────────────────────────────────────
def main():
    excel_path = "Eco_Report_Part1_v1 (1).xlsx"
    df = load_excel(excel_path)
    
    results = run_pipeline(df, n_iter=5000, seed=42)
    if results.empty:
        print("No trees with valid DBH — nothing to plot.")
        return

    # Print and save summary
    summary = results.drop(columns=["AGB_samples"]).round(3)
    print("\n=== AGB summary (Chave et al., 2014; Monte Carlo) ===")
    print(summary.to_string(index=False))
    summary.to_csv("agb_summary.csv", index=False)

    # Generate plots
    make_plots(results, out_prefix="agb")

if __name__ == "__main__":
    main()
