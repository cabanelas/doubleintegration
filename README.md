# Double Integration Analysis

## Pelagic Community Structure: Interannual Variability and Long-Term Change in Pelagic Community Structure Across a Latitudinal Gradient

### Investigating Lagged and Cumulative Biological Responses

This repository contains scripts to test the Double Integration Hypothesis (DIH) proposed by Di Lorenzo & Ohman (2013) across four Long-Term Ecological Research (LTER) sites: California Current Ecosystem (CCE), Northern Gulf of Alaska (NGA), Palmer Station Antarctica (PAL), and Northeast U.S. Shelf (NES). This work is part of the LTER Pelagic Community Structure Working Group, which investigates interannual variability and long-term change in pelagic ecosystems. Specifically, this project examines how marine ecosystems respond to both cyclic and long-term environmental changes using comparative data from multiple LTER sites.

The double integration hypothesis posits that marine populations respond to stochastic environmental forcing through cumulative integration, potentially leading to apparent state changes (Di Lorenzo & Ohman 2013). This analysis aims to identify whether such dynamics are detectable in long-term biological datasets across contrasting ecosystems.

## This repo:
* Perform double integration of normalized biological anomalies.
* Calculate correlations between biological time series and physical drivers.
* Compute autoregressive (AR) coefficients from biological and driver (large-scale indices) time series for the bootstrapping analysis.
* Bootstrap correlations to assess robustness of relationships.

## Run scripts in the following order for best results
Scripts should be run in the following order per site:
1. ARcoef_driver_SITE.R
* Calculates AR(1) coefficients for physical drivers (e.g., PDO, MEI, ONI, AMO).
2. ARcoef_bio_SITE.R
* Calculates AR(1) coefficients for biological time series.
3. doubleintegration_SITE.R
* Performs integrations of normalized biological data.
4. bootstrap_SITE.R
* Bootstrapping.

Replace SITE with the target ecosystem: CCE, PAL, or NGA.
