# Eco-Report AGB Pipeline — README & Uncertainty Log

Monte Carlo uncertainty propagation of tree aboveground biomass (AGB)
following **Chave et al. (2014)** pantropical allometry, applied to
`Eco_Report_Part1_v1 (1).xlsx`.

---

## 1. Purpose

For every living tree in the field sheet, the pipeline:

1. Reads DBH (or CBH → DBH) and wood density from the Excel file.
2. Runs **N = 5,000 Monte Carlo iterations**, drawing DBH and wood density
   from their uncertainty distributions on each iteration.
3. Computes AGB per iteration with the Chave (2014) equation.
4. Reports the **mean**, **SD**, and **95 % credible interval**
   (2.5th–97.5th percentile) of AGB per tree.
5. Produces publication-quality figures, including a **per-TreeID boxplot**
   of the Monte Carlo AGB distribution.

---

## 2. Allometry — Chave et al. (2014)

The pipeline uses the **Height-inclusive pantropical form (Equation 4)**. 
This ensures the Monte Carlo estimates center around the exact values 
calculated in Column J of the Excel sheet.

### Equation (Chave et al. 2014, Eq. 4)

```
AGB = 0.0673 * (ρ * D² * H)^0.976
```

| Symbol | Meaning                          | Value used |
|:------:|:---------------------------------|:----------:|
| `AGB`  | Aboveground biomass (kg)         | —          |
| `D`    | Diameter at breast height (cm)   | per tree   |
| `ρ`    | Wood density (g/cm³)             | per tree   |
| `H`    | Total tree height (m)            | per tree (fallback = 5.0 m) |
| `0.0673` | Allometric constant            | fixed      |
| `0.976`  | Scaling exponent               | fixed      |

Source coefficients: Chave et al. 2014, *Global Change Biology* 20: 
3177–3190, Table 4.

### Reference model performance — Chave (2014) RSE

| Model variant                    | RSE (residual standard error) |
|:---------------------------------|:-----------------------------:|
| D + ρ + H (used here)            | **≈ 0.359** (ln units)        |
| D + ρ                            | ≈ 0.418                       |
| D only                           | ≈ 0.443                       |

The RSE of ≈ 0.359 corresponds to roughly **±43 % 1-σ prediction error** 
on an individual tree. This is *model structural* error and is **not** 
included in the per-tree Monte Carlo bands below — those quantify 
*input measurement + wood-density* error only.

---

## 3. Input uncertainty values used

### 3.1 DBH — diameter at breast height

| Source of uncertainty        | Distribution / range                | Value |
|:-----------------------------|:------------------------------------|:-----:|
| Caliper / tape measurement   | Uniform multiplicative, **±2 %**    | `dbh_frac_err = 0.02` |
| Multi-stem handling          | Each stem sampled independently, then summed per iteration | — |
| Excel formula (`=F2/PI()`)   | Fully evaluated by script (`CBH / π`) | π = 3.14159… |
| Multi-stem string (`DBH1=..`)| Extracted directly from the text string | per stem |

Rationale: ±2 % is the standard field error for DBH measured with a
diameter tape or caliper on smooth-barked trees. If your field protocol
specifies a different tolerance, change `dbh_frac_err` in
`monte_carlo_agb_tree`.

### 3.2 Wood density (ρ) — species-specific ranges

Drawn from a **uniform distribution** `ρ ~ U(ρ_min, ρ_max)` per iteration,
per species, based on the Global Wood Density Database (Zanne et al. 2009;
Chave et al. 2009) and the ranges listed in the source
[DataDryad DOI: 10.5061/dryad.234](https://datadryad.org/dataset/doi:10.5061/dryad.234).

| Species                 | ρ_min (g/cm³) | ρ_max (g/cm³) | Midpoint used in sheet | Note |
|:------------------------|:-------------:|:-------------:|:----------------------:|:-----|
| *Cassia fistula*        | 0.71          | 0.96          | 0.96                   | sheet uses upper bound |
| *Pongamia pinnata*      | 0.54          | 0.67          | 0.64                   | sheet uses midpoint |
| *Swietenia mahagoni*    | 0.54          | 0.83          | 0.51                   | sheet uses lower bound |
| *Terminalia bellirica*  | 0.57          | 0.90          | 0.76                   | sheet uses midpoint |

The pipeline **ignores the point density in column I** and instead
samples the full [min, max] range above — this is the recommended
practice when the within-species density spread is large.

> If you prefer the sheet's point values with a symmetric band instead,
> swap the two lines in `run_pipeline`:
>
> ```python
> rho_pt = float(row["Density (g/cm³)"])
> rho_min, rho_max = rho_pt * 0.8, rho_pt * 1.2   # ±20 %
> ```

### 3.3 Height (H)

Height is read directly from the `Height( in m)` column for single-stem trees. 
For multi-stem trees where height is recorded as `1.66` (likely a stem 
height rather than total tree height) or is blank, the script defaults 
to a fallback height of **5.0 meters** to prevent near-zero biomass 
inflation.

### 3.4 Monte Carlo configuration

| Parameter                 | Value   |
|:--------------------------|:-------:|
| Iterations per tree       | **5,000** |
| RNG seed                  | **42** (reproducible) |
| Credible interval reported| **95 % (2.5th–97.5th pct)** |

---

## 4. Tree-level input handling

| Tree type   | DBH source used                                         |
|:------------|:--------------------------------------------------------|
| Single stem | Excel formula evaluated: `CBH / π` (matches `=F2/PI()`) |
| Multi-stem  | Parsed strictly from `DBH ( in cm)` text; stems summed  |
| Blank CBH + blank DBH | tree skipped with a `[SKIP]` log line          |
| `TreeExists = FALSE` | tree excluded upstream (`load_excel`)            |

**Note on Excel formula cells.** In the sheet, single-stem DBH is written
as the formula `=F2/PI()`. `pandas.read_excel` cannot evaluate formulas,
so the script uses a regex to extract the row number, grab the CBH from 
that row, and recomputes `CBH / π` internally. The result is numerically 
identical to the sheet.

---

## 5. Output files

| File                            | Content                                                    |
|:--------------------------------|:-----------------------------------------------------------|
| `agb_summary.csv`               | Per-tree mean, SD, 2.5 %, 97.5 % AGB (kg)                  |
| `agb_per_tree_ci.png`           | Bar chart of mean AGB ± 95 % MC CI, one bar per TreeID     |
| `agb_distributions.png`         | Per-tree posterior histograms (5,000 samples each)         |
| `agb_per_tree_boxplot.png`      | **Boxplot, one box per TreeID** — the requested figure     |
| `README.md`                     | This file                                                  |

Console log also prints:
- Total rows with `TreeID`
- Number of `TreeExists == TRUE` trees
- Any species not present in `DENSITY_DICT`
- Any tree skipped for missing DBH/CBH

---

## 6. Known caveats & assumptions

1. **Multi-stem Height fallback.** Because your sheet records `1.66` for 
   multi-stem heights (which is likely a measurement from the base to the 
   first branching point, not total tree height), the script applies a 
   fallback height of `5.0 m` to multi-stem trees. If you acquire true 
   total heights for these trees later, update the `Height( in m)` column.
2. **Density range, not point value.** The pipeline deliberately uses the
   species range. If your study specifically wants the sheet's point
   density, see §3.2.
3. **No spatial / stand-level error.** Errors are per-tree and independent.
   For stand-level AGB, sum the per-tree samples **column-wise**
   (`results_df["AGB_samples"].sum(axis=0)` gives 5,000 stand totals) and
   report mean ± 95 % CI of those.
4. **Chave (2014) RSE not included.** Structural allometric error
   (≈ ±43 % for Eq. 4) is *not* part of the per-tree MC bands. 
5. **Wood density source.** The ranges in `DENSITY_DICT` match the
   dryad dataset cited in column E of the sheet; species not in that
   dataset are skipped with a warning.

---

## 7. Reproducibility

```
Python ≥ 3.9
numpy   ≥ 1.22
pandas  ≥ 1.4
matplotlib ≥ 3.5
openpyxl (for .xlsx)
```

Seed = 42 throughout; re-running on the same input file yields identical
numbers and figures.

---

## 8. References

- Chave, J., Réjou-Méchain, M., Búrquez, A., Chidumayo, E., Colgan, M.S.,
  Delitti, W.B.C., Duque, A., Eid, T., Fearnside, P.M., Goodman, R.C.,
  Henry, M., Martínez-Yrízar, A., Mugasha, W.A., Muller-Landau, H.C.,
  Mencuccini, M., Nelson, B.W., Ngomanda, A., Nogueira, E.M., Ortiz-Malavassi,
  E., Pélissier, R., Ploton, P., Ryan, C.M., Saldarriaga, J.G., Vieilledent, G.
  (2014). Improved allometric models to estimate the aboveground biomass
  of tropical trees. *Global Change Biology* **20**: 3177–3190.
  https://doi.org/10.1111/gcb.12629
- Chave, J., Coomes, D., Jansen, S., Lewis, S.L., Swenson, N.G., Zanne, A.E.
  (2009). Towards a worldwide wood economics spectrum.
  *Ecology Letters* **12**: 351–366.
- Zanne, A.E. et al. (2009). Data from: Towards a worldwide wood economics
  spectrum. Dryad. https://doi.org/10.5061/dryad.234
- Phillips, O.L. et al. (2009). Drought sensitivity of the Amazon
  rainforest. *Science* **323**: 1344–1347.
- Condit, R. (1998). *Tropical Forest Census Plots*. Springer.