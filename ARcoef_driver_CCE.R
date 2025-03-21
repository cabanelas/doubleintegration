################################################################################
#############          Pelagic Synthesis           #############################
#############             FEB-2025                 #############################
#############       CCE Double Integration         #############################
## by: Alexandra Cabanelas 
################################################################################
##### CCE LTER
##### Driver = Pacific Decadal Oscillation PDO
##### 1941
## Double Integration Analysis CCE
# Script #1 : ARcoef_driver_CCE

# script to calculate AR coefficient of driver 

## STEP 1 of Monte Carlo analysis
# Step 1 - estimate autoregression coeff from the original data/signals to be
#able to use for creating two red-noise time series 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(forecast) #v8.21; Arima()
library(tseries) #v0.10.54; ADF test
library(astsa) #v2.1; acf2
library(tidyverse)
#library(urca) #v1.3.3

## ------------------------------------------ ##
#            DATA -----
## ------------------------------------------ ##
PDO <- read.csv(file.path("raw", "CCE",
                          "PDO.csv"))

## ------------------------------------------ ##
#            TIDY -----
## ------------------------------------------ ##
# long format
PDO <-  PDO %>%
  pivot_longer(cols = Jan:Dec,
               names_to = "month",
               values_to = "pdo") %>%
  as.data.frame()

# fix date
PDO$Date <- as.Date(paste(PDO$Year, PDO$month, "01"), 
                    format = "%Y %b %d")

PDO$month <- match(PDO$month, month.abb)

PDO <- PDO %>% 
  filter(pdo < 99 & Year > 1940) #filtering 10 yrs before bio data

# create time series object
PDOts <- ts(PDO$pdo, start = c(1940, 1), frequency = 12)

head(PDO)
head(PDOts)

## ------------------------------------------ ##
# 1) Detrend time series (if needed) ----> then AR(1) 
## ------------------------------------------ ##
# A time series that is non-stationary (has a trend or changing variance) 
#should be detrended before AR(1) estimation

## No need to detrend Multivariate ENSO Index TS
## first, visual check...
plot(PDOts, main = "Pacific Decadal Oscillation PDO Index", 
     ylab = "PDO", 
     xlab = "Year", 
     type = "l", 
     col = "black", 
     lwd = 2)
# no obvious upward or downward slope 

PDO_decomp <- decompose(PDOts, type = "additive")  # use "multiplicative" if variance increases over time
plot(PDO_decomp)

# to confirm, can perform Augmented Dickey-Fuller test 
# ADF pval = 0.01 == no need to detrend (if pvalue > 0.5 ts has trend)
adf.test(PDOts) #as we know these tests can be a bit misleading
acf(PDOts)
# double checked anyways and detrending vs not detrending == same AR1 coef

## ------------------------------------------ ##
# 2) Fit AR(1)
## ------------------------------------------ ##
# p = 1 = one autoregressive term
# d = 0 = no differencing 
# q = 0 = no moving average term
# fit AR(1) model
# this is the MLE approach (got same result as OLS)
AR1_model <- Arima(PDOts, order = c(1,0,0))

# print summary
print(summary(AR1_model))
ar1_coefficient <- AR1_model$coef["ar1"]  
ar1_coefficient
## ------------------------------------------ ##
# diagnostic plots
## ------------------------------------------ ##
tsdiag(AR1_model)

# extract residuals
PDO_residuals <- residuals(AR1_model)

# plot original time series with fitted values
par(mfrow = c(2, 2))

ts.plot(PDOts, main = "PDO Time Series with AR(1) Fit", ylab = "PDO")
lines(fitted(AR1_model), col = "red", lwd = 2)

# residuals plot
plot(PDO_residuals, main = "Residuals of AR(1) Model", ylab = "Residuals", type = "l")

# histogram of residuals
hist(PDO_residuals, main = "Histogram of Residuals", xlab = "Residuals")

# Q-Q plot of residuals
qqnorm(PDO_residuals)
qqline(PDO_residuals, col = "red")

# autocorrelation of residuals
par(mfrow = c(2, 1))
acf2(PDO_residuals)  # replaces Acf and pacf

par(mfrow = c(1, 1))

## ------------------------------------------ ##
# 3) Export AR coefficient
## ------------------------------------------ ##
se <- sqrt(diag(vcov(AR1_model)))[1]  # standard error
sigma2 <- AR1_model$sigma2

ar_info_df <- data.frame(ts_name = "PDO",
                         AR_coef = ar1_coefficient,
                         se = se,
                         sigma2 = sigma2)

print(ar_info_df)

#write.csv(ar_info_df, "output/CCE/AR_coef_PDOdriver_CCE.csv")