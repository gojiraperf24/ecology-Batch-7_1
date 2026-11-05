import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Species Wood Density Database [min, max] in g/cm3
DENSITY_DICT = {
    "Cassia fistula": [0.71, 0.96],
    "Pongamia pinnata": [0.54, 0.67],
    "Swietenia mahagoni": [0.54, 0.83],
    "Terminalia bellirica": [0.57, 0.90],
}


def parse_dbh(val):
    if pd.isna(val):
        return []
    if isinstance(val, (int, float)):
        return [float(val)] if val > 0 else []

    numbers = re.findall(r"\d+\.?\d*", str(val))
    return [float(n) for n in numbers if float(n) > 0]


def monty_multi_stem(
    dbh_list, height, rho_min, rho_max, n_draws=10000, seed=42
):
    np.random.seed(seed)

    h_sim = np.random.normal(loc=height, scale=0.10 * height, size=n_draws)
    rho_sim = np.random.uniform(low=rho_min, high=rho_max, size=n_draws)

    chave_rse = 0.357
    agb_total_sim = np.zeros(n_draws)

    for dbh in dbh_list:
        dbh_sim = np.random.normal(loc=dbh, scale=0.3, size=n_draws)
        epsilon = np.random.normal(loc=0, scale=chave_rse, size=n_draws)

        agb_stem = (
            0.0673
            * (rho_sim * (dbh_sim**2) * h_sim) ** 0.976
            * np.exp(epsilon)
        )
        agb_total_sim += agb_stem

    mean_agb = np.mean(agb_total_sim)
    p_2_5 = np.percentile(agb_total_sim, 2.5)
    p_97_5 = np.percentile(agb_total_sim, 97.5)

    return mean_agb, p_2_5, p_97_5


# --- STEP 1: Data Processing Pipeline ---
file_path = r"Eco_Report_Part1_v1 (1).xlsx"
df = pd.read_excel(file_path)

# Rename columns to clean standard names
df.rename(
    columns={
        "Tree ID": "TreeID",  # Rename if your Excel file has a tree ID column
        "DBH ( in cm)": "DBH",
        "Height( in m)": "Height",
        "TreeSpecies": "TreeSpecies",
    },
    inplace=True,
)

# Filter missing essential entries
df.dropna(subset=["TreeSpecies", "Height"], inplace=True)

if "TreeExists" in df.columns:
    df = df[df["TreeExists"] != 0]

results = []

# Enumerate to map rows cleanly to Tree IDs 301 to 325
for i, row in enumerate(df.itertuples()):
    # If Excel already has a TreeID column, use: tree_id = getattr(row, "TreeID", 301 + i)
    tree_id = getattr(row, "TreeID", 301 + i)

    species = row.TreeSpecies
    height = row.Height
    dbh_raw = row.DBH

    if species not in DENSITY_DICT or pd.isna(height) or height <= 0:
        continue

    dbh_stems = parse_dbh(dbh_raw)
    if not dbh_stems:
        continue

    rho_min, rho_max = DENSITY_DICT[species]

    mean_agb, p2_5, p97_5 = monty_multi_stem(
        dbh_list=dbh_stems,
        height=height,
        rho_min=rho_min,
        rho_max=rho_max,
    )

    results.append(
        {
            "Tree_ID": tree_id,  # Inferred tree ID (301 to 325)
            "Species": species,
            "Stem_Count": len(dbh_stems),
            "DBHs_cm": dbh_stems,
            "Height_m": height,
            "Mean_AGB_kg": round(mean_agb, 2),
            "CI_2.5%_kg": round(p2_5, 2),
            "CI_97.5%_kg": round(p97_5, 2),
        }
    )

results_df = pd.DataFrame(results)
print(results_df)

# Save to Excel
results_df.to_excel(
    r"AGB_Results_301_325.xlsx", index=False
)


# --- STEP 2: Plotting Single Figures (Linear & Log Scales) ---
# Sort results by Species and Tree_ID
results_df_sorted = results_df.sort_values(
    by=["Species", "Tree_ID"]
).reset_index(drop=True)

# Calculate asymmetric error bounds
lower_err = (
    results_df_sorted["Mean_AGB_kg"] - results_df_sorted["CI_2.5%_kg"]
)
upper_err = (
    results_df_sorted["CI_97.5%_kg"] - results_df_sorted["Mean_AGB_kg"]
)
yerr = [lower_err, upper_err]

# Custom x-tick labels: "Species Name (Tree ID)"
x_positions = range(len(results_df_sorted))
x_labels = [
    f"{r.Tree_ID}" for r in results_df_sorted.itertuples()
]
n_draws = 10000


def create_agb_plot(yscale="linear", save_path="plot.png"):
    fig, ax = plt.subplots(figsize=(12, 6), dpi=300)

    # Mean point with 2.5th to 97.5th percentile whiskers
    ax.errorbar(
        x=x_positions,
        y=results_df_sorted["Mean_AGB_kg"],
        yerr=yerr,
        fmt="o",
        color="navy",
        ecolor="crimson",
        capsize=4,
        capthick=1.2,
        elinewidth=1.2,
        label="Mean AGB (95% CI)",
    )

    ax.set_yscale(yscale)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=9)

    scale_label = "Log Scale" if yscale == "log" else "Linear Scale"
    ax.set_xlabel("Tree ID", fontsize=11, labelpad=10)
    ax.set_ylabel(
        f"Aboveground Biomass (AGB) [kg] ({scale_label})", fontsize=11
    )
    ax.set_title(
        f"Aboveground Biomass (AGB) Estimates per Tree ({scale_label})",
        fontsize=13,
        pad=12,
    )

    if yscale == "log":
        ax.grid(True, which="both", linestyle="--", alpha=0.4)
    else:
        ax.grid(True, linestyle="--", alpha=0.5)

    ax.legend(loc="upper left")

    # Figure caption describing units and simulation draws count
    caption_text = (
        f"Figure: Aboveground Biomass (AGB) estimates per tree (IDs 301–325) on a {scale_label.lower()}. "
        f"Points denote mean biomass (kg) and whiskers span the 95% confidence interval "
        f"(2.5th to 97.5th percentile) calculated from {n_draws:,} Monte Carlo simulation draws."
    )

    plt.figtext(
        0.5,
        -0.12,
        caption_text,
        wrap=True,
        horizontalalignment="center",
        fontsize=9.5,
        style="italic",
    )

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved successfully to: {save_path}")
    plt.close()


# Save both Linear and Log scale plots
create_agb_plot(
    yscale="linear",
    save_path=r"AGB_Plot_Linear.png",
)

create_agb_plot(
    yscale="log", save_path=r"AGB_Plot_Log.png"
)
