################################################################################
#############          Pelagic Synthesis           #############################
#############             FEB-2025                 #############################
#############       PAL Double Integration         #############################
## by: Alexandra Cabanelas 
################################################################################
##### PAL LTER
##### Driver = Multivariate ENSO Index (MEI)
##### 1983 - 2024
## Double Integration Analysis PAL
# Script #1 : ARcoef_driver_PAL

# script to calculate AR coefficient of driver 

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
MEI <- read.csv(file.path("raw", "PAL",
                              "MEI.csv"))

## ------------------------------------------ ##
#            TIDY -----
## ------------------------------------------ ##
# fix date
MEI$DATE <- as.Date(MEI$DATE, format = "%m/%d/%y")

# extract year and month
MEI$Year <- as.numeric(format(MEI$DATE, "%Y"))
MEI$Month <- as.numeric(format(MEI$DATE, "%m"))

MEI <- MEI %>% 
  filter(Year > 1982) #filtering 10 yrs before bio data 

# create time series object
MEIts <- ts(MEI$MEI, start = c(1983, 1), frequency = 12) 

head(MEI)
head(MEIts)

## ------------------------------------------ ##
# 1) Detrend time series (if needed) ----> then AR(1) 
## ------------------------------------------ ##
# A time series that is non-stationary (has a trend or changing variance) 
#should be detrended before AR(1) estimation

## No need to detrend Multivariate ENSO Index TS
## first, visual check...
plot(MEIts, main = "Multivariate ENSO Index", 
     ylab = "MEI", 
     xlab = "Year", 
     type = "l", 
     col = "black", 
     lwd = 2)
# no obvious upward or downward slope 

MEI_decomp <- decompose(MEIts, type = "additive")  # use "multiplicative" if variance increases over time
plot(MEI_decomp)

# to confirm, can perform Augmented Dickey-Fuller test 
# ADF pval = 0.01 == no need to detrend (if pvalue > 0.5 ts has trend)
adf.test(MEIts) #as we know these tests can be a bit misleading 
acf(MEIts)
# double checked anyways and detrending vs not detrending == same AR1 coef

## ------------------------------------------ ##
# 2) Fit AR(1)
## ------------------------------------------ ##
# p = 1 = one autoregressive term
# d = 0 = no differencing 
# q = 0 = no moving average term
# fit AR(1) model
# this is the MLE approach (got same result as OLS)
AR1_model <- Arima(MEIts, order = c(1,0,0))

# print summary
print(summary(AR1_model))
ar1_coefficient <- AR1_model$coef["ar1"]  
ar1_coefficient
## ------------------------------------------ ##
# diagnostic plots
## ------------------------------------------ ##
tsdiag(AR1_model)

# extract residuals
MEI_residuals <- residuals(AR1_model)

# plot original time series with fitted values
par(mfrow = c(2, 2))

ts.plot(MEIts, main = "MEI Time Series with AR(1) Fit", ylab = "MEI")
lines(fitted(AR1_model), col = "red", lwd = 2)

# residuals plot
plot(MEI_residuals, main = "Residuals of AR(1) Model", ylab = "Residuals", type = "l")

# histogram of residuals
hist(MEI_residuals, main = "Histogram of Residuals", xlab = "Residuals")

# Q-Q plot of residuals
qqnorm(MEI_residuals)
qqline(MEI_residuals, col = "red")

# autocorrelation of residuals
par(mfrow = c(2, 1))
acf2(MEI_residuals)  # replaces Acf and pacf

par(mfrow = c(1, 1))

## ------------------------------------------ ##
# 3) Export AR coefficient
## ------------------------------------------ ##
se <- sqrt(diag(vcov(AR1_model)))[1]  # standard error
sigma2 <- AR1_model$sigma2

ar_info_df <- data.frame(ts_name = "MEI",
                         AR_coef = ar1_coefficient,
                         se = se,
                         sigma2 = sigma2)

print(ar_info_df)

#write.csv(ar_info_df, "output/PAL/AR_coef_MEIdriver_PAL.csv")