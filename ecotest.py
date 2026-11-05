import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ── 1. Wood density dictionary (g/cm³) — min/max per species ───────────────
DENSITY_DICT = {
    "Cassia fistula":        [0.71, 0.96],
    "Pongamia pinnata":      [0.54, 0.67],
    "Swietenia mahagoni":    [0.54, 0.83],
    "Terminalia bellirica":  [0.57, 0.90],
}


# ── 2. Helpers ─────────────────────────────────────────────────────────────
def parse_dbh(val):
    """Return list of positive floats found in a DBH/CBH cell."""
    if val is None:
        return []
    if isinstance(val, (list, tuple, np.ndarray)):
        out = []
        for v in val:
            out.extend(parse_dbh(v))
        return out
    if isinstance(val, (int, float, np.integer, np.floating)):
        if pd.isna(val) or val <= 0:
            return []
        return [float(val)]
    s = str(val).strip()
    if s == "" or s.lower() in ("nan", "na", "none"):
        return []
    return [float(n) for n in re.findall(r"\d+\.?\d*", s) if float(n) > 0]


def chave_agb(dbh_cm, rho_gcm3):
    """Chave et al. (2014) pantropical allometry, D-only form (no H term).
    ln(AGB) = -1.803 - 0.976·ln(ρ) + 3.649·ln(D) - 0.0299·(ln D)²  -> kg
    """
    D   = np.asarray(dbh_cm,  dtype=float)
    rho = np.asarray(rho_gcm3, dtype=float)
    lnD   = np.log(D)
    lnrho = np.log(rho)
    return np.exp(-1.803 - 0.976 * lnrho + 3.649 * lnD - 0.0299 * lnD**2)


def monte_carlo_agb_tree(dbh_cm_list, rho_min, rho_max,
                         n_iter=5000, dbh_frac_err=0.02, seed=42):
    """Per-tree MC: DBH ±2 % uniform, density uniform in [rho_min, rho_max].
    Multi-stem trees: stems summed per iteration."""
    rng = np.random.default_rng(seed)
    dbh = np.asarray(dbh_cm_list, dtype=float)
    if dbh.size == 0:
        return np.zeros(n_iter)
    dbh_s = dbh * rng.uniform(1 - dbh_frac_err, 1 + dbh_frac_err,
                              size=(n_iter, dbh.size))
    rho_s = rng.uniform(rho_min, rho_max, size=(n_iter, dbh.size))
    return chave_agb(dbh_s, rho_s).sum(axis=1)


def tree_dbh_cm(row):
    """DBH (cm) list for one tree row.
    Priority: (1) multi-stem 'DBH1=..' text, (2) numeric DBH,
              (3) CBH / pi (single-stem rows where DBH is an Excel formula)."""
    raw_dbh = row.get("DBH ( in cm)")
    if isinstance(raw_dbh, str) and any(c.isdigit() for c in raw_dbh):
        vals = parse_dbh(raw_dbh)
        if vals:
            return vals
    if isinstance(raw_dbh, (int, float, np.integer, np.floating)) \
            and not pd.isna(raw_dbh) and raw_dbh > 0:
        return [float(raw_dbh)]
    cbh = row.get("CBH(cm)")
    if isinstance(cbh, (int, float, np.integer, np.floating)) \
            and not pd.isna(cbh) and cbh > 0:
        return [float(cbh) / np.pi]
    cbh_vals = parse_dbh(cbh)
    return [v / np.pi for v in cbh_vals]


# ── 3. Excel loader (handles TreeExists as 1.0/0.0 float OR "TRUE"/"FALSE") ─
def load_excel(path):
    df = pd.read_excel(path, sheet_name=0)

    df = df.dropna(how="all")
    df = df[df["TreeID"].notna()].copy()

    def to_bool(v):
        if isinstance(v, (bool, np.bool_)):
            return bool(v)
        if v is None:
            return False
        if isinstance(v, (int, float, np.integer, np.floating)):
            if pd.isna(v):
                return False
            return v != 0
        return str(v).strip().upper() == "TRUE"

    df["TreeExists"] = df["TreeExists"].map(to_bool)
    df["TreeSpecies"] = df["TreeSpecies"].astype(str).str.strip()

    print(f"[load_excel] rows with TreeID  : {len(df)}")
    print(f"[load_excel] TreeExists == True: {int(df['TreeExists'].sum())}")
    missing = sorted(set(df.loc[df['TreeExists'], 'TreeSpecies']) -
                     set(DENSITY_DICT.keys()))
    if missing:
        print(f"[load_excel] species NOT in DENSITY_DICT: {missing}")

    return df[df["TreeExists"]].copy()


# ── 4. Pipeline: one row per tree ─────────────────────────────────────────
def run_pipeline(df, n_iter=5000, seed=42):
    rows = []
    for _, row in df.iterrows():
        species = row["TreeSpecies"]
        if species not in DENSITY_DICT:
            print(f"[SKIP] TreeID {row['TreeID']}: species '{species}' "
                  f"not in DENSITY_DICT")
            continue
        dbh_list = tree_dbh_cm(row)
        if not dbh_list:
            print(f"[SKIP] TreeID {row['TreeID']}: no valid DBH/CBH")
            continue

        rho_min, rho_max = DENSITY_DICT[species]
        samples = monte_carlo_agb_tree(dbh_list, rho_min, rho_max,
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
    ax.set_xticklabels(results_df["TreeID"].astype(int).astype(str),
                       rotation=90)
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
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(3 * n_cols, 2.6 * n_rows),
                             squeeze=False)
    for i, (_, row) in enumerate(results_df.iterrows()):
        ax = axes[i // n_cols, i % n_cols]
        ax.hist(row["AGB_samples"], bins=40, color="#55A868",
                edgecolor="white", alpha=0.85)
        ax.axvline(row["AGB_mean_kg"], color="black", lw=1.2, ls="--",
                   label=f"mean = {row['AGB_mean_kg']:.2f} kg")
        ax.set_title(f"Tree {int(row['TreeID'])} — {row['Species']}",
                     fontsize=8)
        ax.set_xlabel("AGB (kg)", fontsize=8)
        ax.set_ylabel("Count", fontsize=8)
        ax.legend(fontsize=7, loc="upper right")
    for j in range(n, n_rows * n_cols):
        axes[j // n_cols, j % n_cols].set_visible(False)
    fig.suptitle("Monte Carlo AGB distributions per tree", y=1.02,
                 fontsize=11)
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
    bp = ax.boxplot(groups, patch_artist=True, widths=0.6,
                    manage_ticks=False,
                    medianprops={"color": "black"},
                    flierprops={"marker": "o", "markersize": 2, "alpha": 0.4})
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    colors = plt.colormaps["Set3"].resampled(len(groups))(
        np.arange(len(groups)))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_edgecolor("black")
    ax.set_xlabel("Tree ID")
    ax.set_ylabel("Aboveground biomass (kg)")
    ax.set_title("Per-tree AGB uncertainty — Monte Carlo "
                 "(Chave et al., 2014)")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_per_tree_boxplot.png")
    plt.close(fig)

    print(f"Figures saved: {out_prefix}_per_tree_ci.png, "
          f"{out_prefix}_distributions.png, "
          f"{out_prefix}_per_tree_boxplot.png")


# ── 6. Main ────────────────────────────────────────────────────────────────
def main():
    excel_path = "Eco_Report_Part1_v1 (1).xlsx"

    df = load_excel(excel_path)
    print(f"Loaded {len(df)} living tree(s) from {excel_path}")

    results = run_pipeline(df, n_iter=5000, seed=42)
    if results.empty:
        print("No trees with valid DBH — nothing to plot.")
        return

    summary = results.drop(columns=["AGB_samples"]).round(3)
    print("\n=== AGB summary (Chave et al., 2014; Monte Carlo) ===")
    print(summary.to_string(index=False))
    summary.to_csv("agb_summary.csv", index=False)

    make_plots(results, out_prefix="agb")


if __name__ == "__main__":
    main()
