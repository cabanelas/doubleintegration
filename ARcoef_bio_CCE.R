################################################################################
#############          Pelagic Synthesis           #############################
#############             MAR-2025                 #############################
#############        CCE - Double Integration      #############################
## by: Alexandra Cabanelas
################################################################################
## CCE LTER
## Organism = Nyctiphanes simplex
## 1951-2021

# Script #2 : ARcoef_bio_CCE

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
#library(urca) #v1.3.3

## ------------------------------------------ ##
#            Data & Tidy -----
## ------------------------------------------ ##
bio <- read.csv(file.path("raw",
                          "CCE",
                          "nsimplex_CCE.csv")) %>%
  rename(year = Year) %>%
  mutate(
    taxa = "Nsimplex",
    year = as.numeric(year)
  ) %>%
  # --- add the missing years -----
  #no sampling = 1967, 1968, 1971, 1973 and 2020
  #real 0s = 1972, 1976, 2010-2012
  complete(year = seq(min(year), max(year), by = 1)) %>%
  as.data.frame() %>%
  arrange(year) %>%
  mutate(
    #linear interpolation - these analyses dont like NAs
    Abundance = approx(year, Abundance, xout = year)$y,
    # calculate anomalies
    # Abundance = Log10(Abundance per m2 + 1)
    Yc_mean = mean(Abundance, na.rm = TRUE),
    Yc_sd = sqrt(sum((Abundance - Yc_mean)^2, na.rm = TRUE) / sum(!is.na(Abundance))), #pop mean; not sample..
    #z-score
    anomaly_yr = (Abundance - Yc_mean) / Yc_sd
    #could also use scale(Abundance)[,1] if using sample SD
  )

# --- create time series object -----
bioTS <- ts(bio$anomaly_yr, 
            start = min(bio$year), 
            frequency = 1)

## ------------------------------------------ ##
# 1) Detrend time series (if needed) ----> then AR(1)
## ------------------------------------------ ##
# A time series that is non-stationary (has a trend or changing variance) 
#should be detrended before AR(1) estimation
plot(bioTS, main = "Nyctiphanes simplex", 
     ylab = "N. simplex", 
     xlab = "Year", 
     type = "l", 
     col = "black", 
     lwd = 2)

#boxplot(bioTS~cycle(bioTS))
ts.plot(bioTS)
abline(reg=lm(bioTS~time(bioTS)))

acf(bioTS)
# double checked anyways and detrending vs not detrending == almost same AR1 coef

## ------------------------------------------ ##
# 2) Fit AR(1)
## ------------------------------------------ ##
# p = 1 = one autoregressive term
# d = 0 = no differencing 
# q = 0 = no moving average term
# fit AR(1) model
# this is the MLE approach (got same result as OLS)
AR1_model <- Arima(bioTS, order = c(1,0,0))

# print summary
print(summary(AR1_model))
ar1_coefficient <- AR1_model$coef["ar1"]  
ar1_coefficient

#sarima(bioTS, 1, 0, 0) #same thing, but gives various diag plots

## ------------------------------------------ ##
# diagnostic plots
## ------------------------------------------ ##
tsdiag(AR1_model) #forecast::checkresiduals(AR1_model)

# extract residuals
bio_residuals <- residuals(AR1_model)

# plot original time series with fitted values
par(mfrow = c(2, 2))

ts.plot(bioTS, main = "N. simplex Time Series with AR(1) Fit", ylab = "N. simplex")
lines(fitted(AR1_model), col = "red", lwd = 2)

# residuals plot
plot(bio_residuals, main = "Residuals of AR(1) Model", ylab = "Residuals", type = "l")

# histogram of residuals
hist(bio_residuals, main = "Histogram of Residuals", xlab = "Residuals")

# Q-Q plot of residuals
qqnorm(bio_residuals)
qqline(bio_residuals, col = "red")

# autocorrelation of residuals
par(mfrow = c(2, 1))
acf2(bio_residuals)  # replaces Acf and pacf
# residuals behave like white noise 
par(mfrow = c(1, 1))

## ------------------------------------------ ##
# 3) Export AR coefficient
## ------------------------------------------ ##
se <- sqrt(diag(vcov(AR1_model)))[1]  # standard error
sigma2 <- AR1_model$sigma2

ar_info_df <- data.frame(site = "CCE",
                         spp = "N_simplex",
                         AR_coef = ar1_coefficient,
                         se = se,
                         sigma2 = sigma2)

print(ar_info_df)
#write.csv(ar_info_df, "output/CCE/ARcoef_Nsimplexbio_CCE.csv")