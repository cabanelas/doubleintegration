################################################################################
#############        LTER Pelagic Synthesis WG     #############################
#############        NGA - Double Integration      #############################
#############      Biology AR(1) coefficient       #############################
## by: Alexandra Cabanelas
## created MAR-2025, updated OCT-2025
################################################################################
## Northern Gulf of Alaska LTER
## Organism = Neocalanus cristatus = SPRING
## 1998 - 2022

# Script #2 : 02_ARcoef_NcrisSprbio_NGA

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
#            Data & Tidy -----
## ------------------------------------------ ##
bio <- read.csv(file.path("raw",
                          "NGA",
                          "NeocalanusCristBiomass_Sprv2.csv")) %>%
  rename(year = Year) %>%
  mutate(
    # ---remove S from year col -----
    year = as.numeric(gsub("S", "", year)),
    # --- calculate anomalies -----
    Yc_mean = mean(LogMean, na.rm = TRUE),
    Yc_sd = sd(LogMean, na.rm = TRUE),
    anomaly_yr = (LogMean - Yc_mean) / Yc_sd
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
plot(bioTS, main = "Neocalanus Cristatus - spring", ylab = "N. cristatus", 
     xlab = "Year", type = "l", col = "black", lwd = 2)

#boxplot(bioTS~cycle(bioTS))
ts.plot(bioTS); abline(reg=lm(bioTS~time(bioTS)))

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
#sarima(bioTS, 1, 0, 0) #same thing, but gives various diag plots

# print summary
print(summary(AR1_model))
ar1_coefficient <- AR1_model$coef["ar1"]  
ar1_coefficient #0.3562206 

## ------------------------------------------ ##
# diagnostic plots
## ------------------------------------------ ##
tsdiag(AR1_model) #forecast::checkresiduals(AR1_model)

# extract residuals
bio_residuals <- residuals(AR1_model)

# plot original time series with fitted values
par(mfrow = c(2, 2))

ts.plot(bioTS, main = "N. cristatus Time Series with AR(1) Fit", ylab = "N. cristatus")
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
par(mfrow = c(1, 1))
#these dont look great?

## ------------------------------------------ ##
# 3) Export AR coefficient
## ------------------------------------------ ##
se <- sqrt(diag(vcov(AR1_model)))[1]  # standard error
sigma2 <- AR1_model$sigma2
n <- nrow(bio)

ar_info_df <- data.frame(site = "NGA",
                         spp = "NeocalanusCrist",
                         season = "spring",
                         AR_coef = ar1_coefficient,
                         se = se,
                         sigma2 = sigma2,
                         n = n)

print(ar_info_df)
#write.csv(ar_info_df, "output/NGA/ARcoef_NeoCriSpbio_NGA.csv")