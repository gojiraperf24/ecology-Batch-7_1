# Eco-Report AGB Pipeline — README

This document explains the Python code we wrote for our B.Tech Semester 1 Ecology project. The code calculates the Aboveground Biomass (AGB) of trees using the Chave et al. (2014) equation and figures out the error margins using a Monte Carlo simulation.

## 1. What Does The Code Do?
For every living tree in our Excel sheet (`Eco_Report_Part1_v1 (1).xlsx`), the script does the following:
1. Reads the tree's diameter (DBH) and wood density.
2. Runs a simulation **5,000 times** per tree, changing the DBH, height, and density slightly each time to see how the math changes.
3. Calculates the biomass for every single run.
4. Finds the average, standard deviation, and the 95% confidence range (the 2.5th and 97.5th percentiles) for each tree.
5. Makes graphs (including a cool boxplot) and saves an Excel file with all the final stats.

## 2. The Formula We Used
We used the Chave (2014) Equation 4, which includes Height:
**AGB = 0.0673 * (Density * DBH² * Height)^0.976**

### About the Model Error
According to Chave's paper, this equation isn't perfect. It has a "Residual Standard Error" (RSE) of about **0.359**. This basically means that even if we measure a tree perfectly, the formula itself has a natural error of about **±43%** for a single tree. 
*Note: We added this 43% error into our simulation using a lognormal curve, so the final error margins are very realistic.*

## 3. How We Handled the Inputs and Errors
Instead of just putting one number into the formula, we told the code to guess a slightly different number 5,000 times based on real-world error rates.

**3.1 DBH (Diameter)**
*   We assumed a measurement error of **±2.5%**. The code picks a random number from a normal distribution around our measured DBH.
*   For trees with multiple stems (like "DBH1=0.767 DBH2=0.318"), the code calculates the biomass of each stem separately and adds them together.
*   Some DBH cells in our Excel sheet were formulas (like `=F2/PI()`). Since Python can't read Excel formulas, we wrote a small piece of code to just calculate `CBH / Pi` for us.

**3.2 Wood Density**
*   Instead of just using one flat number, we told the code to use a **15% coefficient of variance (CV)** and pick a random density from a lognormal distribution.
*   We used these specific average densities for our species:
    *   *Cassia fistula*: 0.96
    *   *Pongamia pinnata*: 0.64
    *   *Swietenia mahagoni*: 0.51
    *   *Terminalia bellirica*: 0.76

**3.3 Height**
*   We gave Height a **10% CV** using a normal distribution.
*   **Important:** If a tree had a blank or invalid height, we just skipped that tree entirely. We didn't want to guess a fake height.

**3.4 Simulation Settings**
*   Iterations: **5,000 per tree**.
*   Seed: **42** (this just makes sure the random numbers are the same every time we run the code).
*   Confidence Interval: **95%**.

## 4. What Happens to Missing or Bad Data?
*   If `TreeExists` is FALSE, the tree is ignored.
*   If the tree is a single stem, we calculate DBH using the CBH formula.
*   If the tree has a blank height or blank DBH, the code prints a `[SKIP]` message in the terminal and moves on to the next tree.

## 5. What Files Does It Make?
When you run the code, it creates these files:
*   `agb_summary.csv`: A simple table with the average, SD, and 95% range for each tree.
*   `agb_detailed_results.xlsx`: A detailed Excel file with the mean, median, 2.5th percentile, 97.5th percentile, and the CV% for each tree.
*   `agb_per_tree_ci.png`: A bar graph showing the average biomass and error bars.
*   `agb_distributions.png`: A bunch of little histograms showing the 5,000 simulation results for each tree.
*   `agb_per_tree_boxplot.png`: A clean boxplot (we hid the extreme outlier dots so it looks neat and readable).

## 6. Things to Keep in Mind (Caveats)
1.  **No fake heights:** Because we skipped trees with blank heights, our final tree count might be slightly lower than the total in the Excel sheet.
2.  **Stand-level totals:** The errors are calculated per tree. If we want to find the error for the whole forest, we just add up the 5,000 simulation results for every tree column by column.
3.  **Chave Error is included:** Unlike some basic formulas, we actually included the ±43% model error in our lognormal simulation, which is why our CV% is around 42%.

## 7. How to Run It
You just need Python installed with these basic libraries:
*   numpy
*   pandas
*   matplotlib
*   openpyxl (for reading/writing Excel files)

Just change the file path in the `main()` function to point to our Excel file, and run it!