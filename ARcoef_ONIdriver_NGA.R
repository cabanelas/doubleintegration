################################################################################
#############          Pelagic Synthesis           #############################
#############             FEB-2025                 #############################
#############       NGA Double Integration         #############################
## by: Alexandra Cabanelas 
################################################################################
##### NGA LTER
##### Driver = Oceanic Nino Index = ONI
##### 1988
## Double Integration Analysis NGA
# Script #1 : ARcoef_ONIdriver_NGA

# script to calculate AR coefficients of driver 

## STEP 1 of Monte Carlo analysis
# Step 1 - estimate autoregression coeff from the original data/signals to be
#able to use for creating two red-noise time series 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(forecast) #v8.21; Arima()
library(tseries) #v0.10.54; ADF test
library(astsa) #v2.1; acf2 (optional)
library(tidyverse)
#library(urca) #v1.3.3

## ------------------------------------------ ##
#            DATA -----
## ------------------------------------------ ##
ONI <- read.csv(file.path("raw", "NGA",
                          "ONI.csv"))

## ------------------------------------------ ##
#            TIDY -----
## ------------------------------------------ ##
ONI <- ONI %>%
  mutate(month = case_when(
    SEAS == "DJF" ~ 1,
    SEAS == "JFM" ~ 2,
    SEAS == "FMA" ~ 3,
    SEAS == "MAM" ~ 4,
    SEAS == "AMJ" ~ 5,
    SEAS == "MJJ" ~ 6,
    SEAS == "JJA" ~ 7,
    SEAS == "JAS" ~ 8,
    SEAS == "ASO" ~ 9,
    SEAS == "SON" ~ 10,
    SEAS == "OND" ~ 11,
    SEAS == "NDJ" ~ 12
  ))

# fix date
#ONI$YR <- as.numeric(ONI$YR)
ONI$Date <- as.Date(paste(ONI$YR, ONI$month, "01"), 
                    format = "%Y %m %d")

ONI <- ONI %>% 
  filter(YR > 1987) #filtering 10 yrs before bio data 

# create time series object
ONIts <- ts(ONI$ANOM, start = c(1988, 1), frequency = 12)

head(ONI)
head(ONIts)

## ------------------------------------------ ##
# 1) Detrend time series (if needed) ----> then AR(1) 
## ------------------------------------------ ##
# A time series that is non-stationary (has a trend or changing variance) 
#should be detrended before AR(1) estimation

## No need to detrend Multivariate ENSO Index TS
## first, visual check...
plot(ONIts, main = "ONI Index", 
     ylab = "ONI", 
     xlab = "Year", 
     type = "l", 
     col = "black", 
     lwd = 2)
# no obvious upward or downward slope 

ONI_decomp <- decompose(ONIts, type = "additive")  # use "multiplicative" if variance increases over time
plot(ONI_decomp)

# to confirm, can perform Augmented Dickey-Fuller test 
# ADF pval = 0.01 == no need to detrend (if pvalue > 0.5 ts has trend)
adf.test(ONIts) #as we know these tests can be a bit misleading
acf(ONIts)
# double checked anyways and detrending vs not detrending == same AR1 coef

## ------------------------------------------ ##
# 2) Fit AR(1)
## ------------------------------------------ ##
# p = 1 = one autoregressive term
# d = 0 = no differencing 
# q = 0 = no moving average term
# fit AR(1) model
# this is the MLE approach (got same result as OLS)
AR1_model <- Arima(ONIts, order = c(1,0,0))

# print summary
print(summary(AR1_model))
ar1_coefficient <- AR1_model$coef["ar1"]  
ar1_coefficient
## ------------------------------------------ ##
# diagnostic plots
## ------------------------------------------ ##
tsdiag(AR1_model)

# extract residuals
ONI_residuals <- residuals(AR1_model)

# plot original time series with fitted values
par(mfrow = c(2, 2))

ts.plot(ONIts, main = "ONI Time Series with AR(1) Fit", ylab = "ONI")
lines(fitted(AR1_model), col = "red", lwd = 2)

# residuals plot
plot(ONI_residuals, main = "Residuals of AR(1) Model", ylab = "Residuals", type = "l")

# histogram of residuals
hist(ONI_residuals, main = "Histogram of Residuals", xlab = "Residuals")

# Q-Q plot of residuals
qqnorm(ONI_residuals)
qqline(ONI_residuals, col = "red")

# autocorrelation of residuals
par(mfrow = c(2, 1))
acf2(ONI_residuals)  # replaces Acf and pacf

par(mfrow = c(1, 1))

## ------------------------------------------ ##
# 3) Export AR coefficient
## ------------------------------------------ ##
se <- sqrt(diag(vcov(AR1_model)))[1]  # standard error
sigma2 <- AR1_model$sigma2

ar_info_df <- data.frame(ts_name = "ONI",
                         AR_coef = ar1_coefficient,
                         se = se,
                         sigma2 = sigma2)

print(ar_info_df)

#write.csv(ar_info_df, "output/NGA/AR_coef_ONIdriver_NGA.csv")