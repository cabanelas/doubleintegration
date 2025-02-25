# doubleintegration

# Double Integration Analysis

This repository contains code and analyses for testing the Double Integration Hypothesis (DIH) proposed by Di Lorenzo & Ohman (2013) across four Long-Term Ecological Research (LTER) sites: California Current Ecosystem (CCE), Northern Gulf of Alaska (NGA), Palmer Station Antarctica (PAL), and Northeast U.S. Shelf (NES). This work is part of the LTER Pelagic Community Structure Working Group, which investigates interannual variability and long-term change in pelagic ecosystems. Specifically, this project examines how marine ecosystems respond to both cyclic and long-term environmental changes using comparative data from multiple LTER sites.



## Scripts 

1.	Anomalies_bio  
a.	To calculate zscore anomalies  
b.	Output: zscore_anomalies_3taxazp_28MAY  

2.	calculateARcoefficient_drivers_v3  
a.	detrending time series   
b.	calculating AR coefficient   
c.	Output: detrended_driver_time_series.csv; AR_coef_drivers_allAR1_afterDetrend.csv  

3.	calculateARcoefficient_bio_v3  
a.	detrending time series   
b.	calculating AR coefficient   
c.	Output: detrended_BIOLOGY_time_series.csv; AR_coef_bio_allAR1_afterDetrend_wSampleSize.csv  

4.	Integration_Corr_customTAU_DetrendedBioDriver_v2  
a.	Integration, correlation, pvalues, plots  
b.	Output: IntCorrelations_DetrendBioDriver_6090240TAU.csv  
c.	Figures  

5.	Boot_parallel_customTAU  
a.	Bootstrapping  
b.	Output: results_boot_6090240TAU_RAW.csv  

6.	Final_plots_Conv_Boot_pvalues  
a.	Final plots with conventional and bootstrapped pvalues  
