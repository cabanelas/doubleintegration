################################################################################
#############          Pelagic Synthesis           #############################
#############             FEB-2025                 #############################
#############       CCE - Double Integration       #############################
## by: Alexandra Cabanelas 
################################################################################
## CCE LTER
## Driver = Pacific Decadal Oscillation PDO
## 1941

# Script #1 : ARcoef_PDOdriver_CCE

# script to calculate AR coefficient of driver 

## STEP 1 of Monte Carlo analysis
# Step 1 - estimate autoregression coeff from the original data/signals to be
#able to use for creating two red-noise time series 
# the goal is to use the AR1 coefficient to create red-noise surrogates that mimic
# the autocorreltion structure of the original PDO 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(forecast) #v8.21; Arima()
library(tseries) #v0.10.54; ADF test
library(astsa) #v2.1; acf2 (optional)
library(tidyverse)
#library(urca) #v1.3.3

## ------------------------------------------ ##
#            Data & Tidy -----
## ------------------------------------------ ##
PDO <- read.csv(file.path("raw", 
                          "CCE",
                          "PDO.csv")) %>% # contains up to Dec 2024
       # long format
       pivot_longer(cols = Jan:Dec, 
                    names_to = "month", 
                    values_to = "pdo") %>%
       as.data.frame() %>%
       mutate(
         # fix date
         Date = as.Date(paste(Year, month, "01"),
                        format = "%Y %b %d"),
         #turn month into integer
         month = match(month, month.abb)
       ) %>%
       arrange(Date) %>%
       #filtering 10 yrs before bio data
       filter(Year >= 1941, Year <= 2021, pdo < 99) 

# --- create time series object -----
PDOts <- ts(PDO$pdo, start = c(1941, 1), frequency = 12)

head(PDO)
head(PDOts)

## ------------------------------------------ ##
# 1) Detrend time series (if needed) ----> then AR(1) 
## ------------------------------------------ ##
# A time series that is non-stationary (has a trend or changing variance) 
#should be detrended before AR(1) estimation

## No need to detrend PDO Index TS
## visual check...
plot(PDOts, main = "Pacific Decadal Oscillation PDO Index", 
     ylab = "PDO", 
     xlab = "Year", 
     type = "l", 
     col = "black", 
     lwd = 2)
# no obvious upward or downward slope 
# PDO is already an anomaly index so unlikely that detrending is needed

PDO_decomp <- decompose(PDOts, type = "additive")  # use "multiplicative" if variance increases over time
plot(PDO_decomp)

# to confirm, can perform Augmented Dickey-Fuller test 
# ADF pval = 0.01 == no need to detrend
#if pvalue > 0.05 ts has trend; fail to reject unit root
adf.test(PDOts) #these tests can be a bit misleading
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
summary(AR1_model)
ar1_coefficient <- AR1_model$coef["ar1"] #or coef(AR1_model)[["ar1"]]
ar1_coefficient
## ------------------------------------------ ##
# diagnostic plots
## ------------------------------------------ ##
tsdiag(AR1_model) #forecast::checkresiduals(AR1_model)

# extract residuals
PDO_residuals <- residuals(AR1_model)

# plot original time series with fitted values
par(mfrow = c(2, 2))

# ts (top left)
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
se <- sqrt(diag(vcov(AR1_model)))[1]  # standard error. or  sqrt(vcov(AR1_model)["ar1","ar1"])
sigma2 <- AR1_model$sigma2
# sd <- sqrt(sigma2)
# ci <- ar1_coefficient + c(-1.96, 1.96) * se
ar_info_df <- data.frame(ts_name = "PDO",
                         AR_coef = ar1_coefficient,
                         se = se,
                         sigma2 = sigma2,
                         n = n) #time series length

print(ar_info_df)

#write.csv(ar_info_df, "output/CCE/ARcoef_PDOdriver_CCE.csv")