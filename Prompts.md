\# AI Prompts Used to Generate the Python Script



I used the following prompts to generate the file ecotest.py. Using the given stepped prompts, the AI was able to build the data pipeline, the Monte Carlo simulation, and the plotting functions accurately.



\### Data Loading and Cleaning



```text

I need a Python script to calculate tree aboveground biomass using the Chave 2014 equation. I have an Excel file with tree data. Write the code to load it and clean it up. 

It needs to ignore trees where tree exists is FALSE. For multi-stem trees, the DBH is written as text like "DBH1=0.767 DBH2=0.318", so it needs to parse those out and sum the stems. if a tree's height is blank or invalid, just skip that tree completely.

Add a Monte Carlo simulation to the pipeline (5,000 runs per tree) to calculate error margins. Use Normal distribution for DBH with CV of 2.5%, Normal for Height with 15% CV, Wood density with lognormal distribution and 15% CV and Chave model error of sigma 0.359 with lognormal distribution.


For multi-stem trees simulate each stem separately and add them together for each run


Lastly, write the matplotlib code to save three PNG graphs:

1\\. A bar chart of the mean biomass per tree with 95% error bars.

2\\. A grid of histograms showing the 5,000 simulated samples for each tree.

3\\. A boxplot of the distributions. But please hide the outliers in the boxplot (use showfliers=False) so we don't see the crazy tail dots from the lognormal distribution, just make the boxes look clean.


