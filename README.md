\# Eco-Report AGB Pipeline — README \& Uncertainty Log



Monte Carlo uncertainty propagation of tree aboveground biomass (AGB)

following \*\*Chave et al. (2014)\*\* pantropical allometry, applied to

`Eco\_Report\_Part1\_v1 (1).xlsx`.



\---



\## 1. Purpose



For every living tree in the field sheet, the pipeline:



1\. Reads DBH (or CBH → DBH) and wood density from the Excel file.

2\. Runs \*\*N = 5,000 Monte Carlo iterations\*\*, drawing DBH and wood density

&#x20;  from their uncertainty distributions on each iteration.

3\. Computes AGB per iteration with the Chave (2014) equation.

4\. Reports the \*\*mean\*\*, \*\*SD\*\*, and \*\*95 % credible interval\*\*

&#x20;  (2.5th–97.5th percentile) of AGB per tree.

5\. Produces publication-quality figures, including a \*\*per-TreeID boxplot\*\*

&#x20;  of the Monte Carlo AGB distribution.



\---



\## 2. Allometry — Chave et al. (2014)



The pipeline uses the \*\*D-only (no height) pantropical form\*\* because

tree height in the sheet is recorded only for single-stem trees and is

missing/unreliable for multi-stem individuals.



\### Equation (Chave et al. 2014, Eq. 7)



