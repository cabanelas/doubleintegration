################################################################################
#############          Pelagic Synthesis           #############################
#############             MAR-2025                 #############################
#############        CCE -  Double Integration     #############################
## by: Alexandra Cabanelas
################################################################################
## Double Integration Analysis CCE
# Script #3 : ARcoef_bio_CCE

# script to calculate AR coefficient of the biology time series
# N. simplex 

## STEP 1 of Monte Carlo analysis
# Step 1 - estimate autoregression coeff from the original data/signals to be
#able to use for creating two red-noise time series 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(tidyverse) #v2.0.0
library(here) #v1.0.1
library(forecast) #v8.21
library(tseries) #v0.10.54; ADF test
library(urca) #v1.3.3
#library(tibble) #v3.2.1
#library(purrr) #v1.0.2
#library(tidyr) #v1.3.1
library(astsa) #v2.1

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
bioTS <- read.csv(file.path("raw",
                            "CCE",
                            "Euphausiids_CCE.csv")) %>%
  rename(year = Year,
         anomaly_yr = Anomaly_yr)

## ------------------------------------------ ##
#            Tidy -----
## ------------------------------------------ ##
bioTS$date <- as.Date(paste0(bioTS$year, "-03-01"))
#bioTS$year <- format(bioTS$date, "%Y")
bioTS$year <- as.numeric(bioTS$year)

#need to add the missing years 
#no sampling 1967, 1968, 1971, 1973 and 2020
#real 0s = 1972, 1976, 2010-2012
full_years <- data.frame(year = seq(min(bioTS$year), max(bioTS$year), by = 1))
bioTS <- full_years %>%
  left_join(bioTS, by = "year")

#linear interpolation 
#these analyses dont like NAs
bioTS$anomaly_yr <- approx(bioTS$year, 
                           bioTS$anomaly_yr,
                           xout = bioTS$year)$y

## ------------------------------------------ ##
#            Make ts object -----
## ------------------------------------------ ##
bioTS_ts <- ts(bioTS$anomaly_yr, 
               start = min(bioTS$year), 
               frequency = 1)

## ------------------------------------------ ##
# 1) Detrend time series (if needed) ----> then AR(1)
## ------------------------------------------ ##
# A time series that is non-stationary (has a trend or changing variance) 
#should be detrended before AR(1) estimation
plot(bioTS_ts, main = "Limacina", 
     ylab = "Limacina", 
     xlab = "Year", 
     type = "l", 
     col = "black", 
     lwd = 2)

boxplot(bioTS_ts~cycle(bioTS_ts))

ts.plot(bioTS_ts)
abline(reg=lm(bioTS_ts~time(bioTS_ts)))

acf(bioTS_ts)
# double checked anyways and detrending vs not detrending == almost same AR1 coef

## ------------------------------------------ ##
# 2) Fit AR(1)
## ------------------------------------------ ##
# p = 1 = one autoregressive term
# d = 0 = no differencing 
# q = 0 = no moving average term
# fit AR(1) model
# this is the MLE approach (got same result as OLS)
AR1_model <- Arima(bioTS_ts, order = c(1,0,0))

# print summary
print(summary(AR1_model))
ar1_coefficient <- AR1_model$coef["ar1"]  
ar1_coefficient

#sarima(bioTS_ts, 1, 0, 0) #same thing, but gives various diag plost
## ------------------------------------------ ##
# diagnostic plots
## ------------------------------------------ ##
tsdiag(AR1_model) 

# extract residuals
bio_residuals <- residuals(AR1_model)

# plot original time series with fitted values
par(mfrow = c(2, 2))

ts.plot(bioTS_ts, main = "N. simplex Time Series with AR(1) Fit", ylab = "N. simplex")
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
#write.csv(ar_info_df, "output/CCE/AR_coef_bio_CCE.csv")
