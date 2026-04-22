Copula-Based Bias Correction for Wind Speed Data
This repository provides code to reproduce the methods and results from the paper:
A Copula-Based Bivariate Bias Correction Method for Super-Resolved Physics-Based Wind Speed Data
Published in IEEE
Overview
Accurate wind speed data is essential for applications such as renewable energy integration, power system resilience, and extreme event analysis. Physics-based and reanalysis datasets (e.g., ERA5 and MERRA-2) often exhibit systematic biases relative to observations, particularly in the tails of the distribution.
This repository implements a copula-based bivariate bias correction framework designed to:
- Correct marginal distribution biases in wind speed data
- Preserve dependence structure between variables
- Improve representation of extreme wind events
- Enable consistent integration across multiple datasets
Unlike univariate bias correction methods, this approach explicitly models joint dependence using copulas, leading to more physically consistent corrected data.
