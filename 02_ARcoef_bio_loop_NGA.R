################################################################################
#############        LTER Pelagic Synthesis WG     #############################
#############        NGA - Double Integration      #############################
#############      Biology AR(1) coefficient Loop  #############################
## by: Alexandra Cabanelas
## created OCT-2025
################################################################################
## Northern Gulf of Alaska LTER
## Looping all taxa

# Script #2 : 02_ARcoef_bio_loop_NGA

# script to calculate AR coefficient of the biology time series

## STEP 1 of Monte Carlo analysis
# Step 1 - estimate autoregression coeff from the original data/signals to be
#able to use for creating two red-noise time series 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(tidyverse) #v2.0.0
library(here) #v1.0.1
library(forecast) #v8.21; Arima()
library(tseries) #v0.10.54; ADF test
library(astsa) #v2.1; acf2 (optional)

## ------------------------------------------ ##
#           Data -----
## ------------------------------------------ ##
bio <- read_csv(file.path("raw",
                          "NGA",
                          "allzooplankton_NGA.csv")) %>% #created in NGA_zooplankton.R
  select(-1)

## ------------------------------------------ ##
# Calculate AR coefficient for all bio/taxa data -----
## ------------------------------------------ ##
ar_info_df <- bio %>%
  group_by(taxa, season) %>%
  group_split() %>%
  map_df(function(df_group) {
    # ensure data is sorted by year
    df_group <- df_group %>% arrange(Year)
    
    # extract metadata
    this_taxa   <- unique(df_group$taxa)
    this_season <- unique(df_group$season)
    year_range  <- range(df_group$Year)
    n_years     <- nrow(df_group)
    
    # --- create time series object -----
    bioTS <- ts(df_group$Anomaly_yr,
                start = min(df_group$Year),
                frequency = 1)
    
    # --- fit AR(1) model -----
    AR1_model <- Arima(bioTS, order = c(1, 0, 0)) #MLE approach
    
    # --- extract info -----
    ar1_coefficient <- AR1_model$coef["ar1"]
    se              <- sqrt(diag(vcov(AR1_model)))[1]
    sigma2          <- AR1_model$sigma2
    
    tibble(
      site       = "NGA",
      taxa       = this_taxa,
      season     = this_season,
      year_start = year_range[1],
      year_end   = year_range[2],
      n          = n_years,
      AR_coef    = ar1_coefficient,
      se         = se,
      sigma2     = sigma2
    )
  })

print(ar_info_df)

#write.csv(ar_info_df, "output/NGA/ARcoef_All_bio_NGA.csv")