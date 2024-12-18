################################################################################
#############          Pelagic Synthesis           #############################
#############             MAR-2024                 #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas 
################################################################################
## Double Integration Analysis NES
# Script #2 : calculateARcoefficient_drivers

# script to detrend time series and calculate AR coefficients of drivers 

## STEP 1 of Monte Carlo analysis
# Step 1 - estimate autoregression coeff from the original data/signals to be
#able to use for creating two red-noise time series 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(forecast) #v8.21
library(tseries) #v0.10.54; ADF test
library(urca) #v1.3.3
library(astsa) #v2.1

## ------------------------------------------ ##
#            DATA -----
## ------------------------------------------ ##

# List of data frame names 
data_frames <- c("AMO", "NAO", "AO", "GSI")
ts_names <- c("AMOts", "NAOts", "AOts", "GSIts")

# Create empty lists
data_frame_list <- list()
ts_list <- list()

# Loop through each data frame and create corresponding time series object
for (i in seq_along(data_frames)) {
  # Read data - CSV files
  df <- read.csv(paste("output/", data_frames[i], "Integrations.csv", sep = ""), 
                 header = TRUE)
  
  # Tidy data
  df <- df[, c(2:4)]
  df$Year <- as.numeric(substr(df$time, 1, 4))
  df$month <- factor(substr(df$time, 6, 7), levels = sprintf("%02d", 1:12))
  
  data_frame_list[[data_frames[i]]] <- df
  
  # Create time series object
  assign(ts_names[i], ts(df[, 2], start = c(1973, 1), frequency = 12)) #using raw driver not normalized
  ts_list[[ts_names[i]]] <- ts(df[, 2], start = c(1973, 1), frequency = 12)
  assign(data_frames[i], df)
}
#change to ts(df[, 3] for normalized values of drivers

# List of time series and their names
ts_list <- list(AMO = AMOts, NAO = NAOts, AO = AOts, GSI = GSIts)

ts_names <- c("AMO", "NAO", "AO", "GSI")
#AMOts <- ts_list[["AMO"]]
#NAOts <- ts_list[["NAO"]]
#AOts <- ts_list[["AO"]]
#GSIts <- ts_list[["GSI"]]


## ------------------------------------------ ##
# 1) Detrend time series ----> then AR(1) 
## ------------------------------------------ ##
# not all my time series are stationary so to be able to apply AR(1) model, I detrended them all

#de-trend all the time series by subtracting the simple linear regression line. 
#That way you don't have to think about stationarity at all. 
#All time series will be stationary and you preserve all the fluctuations 
#you're interested in

# decomposition -> reduces time series into 3 components
# trend, seasonal effects, random errors 
# model the random errors as stationary process 
#an alternative to decomposition is differencing - for removing trends 

# Function to detrend a time series using linear regression
detrend_time_series <- function(ts) {
  time <- 1:length(ts)
  lm_trend <- lm(ts ~ time) # Fit linear model
  fitted_trend <- lm_trend$fitted.values # Extract fitted values (trend component)
  #fitted_values <- predict(lm_trend)
  detrended_ts <- ts - fitted_trend
  residuals <- lm_trend$residuals
  return(list(detrended_ts = detrended_ts, 
              residuals = residuals,
              fitted_trend = fitted_trend))
}

#detrended_ts_amo <- detrend_time_series(ts_list[[1]])
#detrended_ts_nao <- detrend_time_series(ts_list[[2]])
#detrended_ts_ao <- detrend_time_series(ts_list[[3]])
#detrended_ts_gsi <- detrend_time_series(ts_list[[4]])

detrended_results_list <- lapply(ts_list, detrend_time_series)

detrended_ts_list <- lapply(detrended_results_list, function(x) x$detrended_ts)
residuals_list <- lapply(detrended_results_list, function(x) x$residuals)
fit_trend_list <- lapply(detrended_results_list, function(x) x$fitted_trend)

# save detrended ts
AMO_detrended <- detrended_ts_list$AMO
NAO_detrended <- detrended_ts_list$NAO
AO_detrended <- detrended_ts_list$AO
GSI_detrended <- detrended_ts_list$GSI

merged_df <- cbind(AMO_detrended, NAO_detrended, AO_detrended, GSI_detrended)
merged_df_1 <- as.data.frame(merged_df)
start_date <- as.Date("1973-01-15")
num_months <- nrow(merged_df)  
dates <- seq(start_date, by = "month", length.out = num_months)

# Add the time column to merged_df_1
merged_df_1 <- cbind(time = dates, merged_df_1)


#write.csv(merged_df_1, "detrended_driver_time_series.csv", row.names = TRUE)

## ------------------------------------------ ##
# PLOTS
## ------------------------------------------ ##

par(mfrow = c(2,1))
for (i in seq_along(ts_list)) {
  original_ts <- ts_list[[i]]
  detrended_ts <- detrended_ts_list[[i]]
  
  # Plot original time series
  plot(original_ts, main = paste("Original", names(ts_list)[i]), 
       ylab = names(ts_list)[i], xlab = "Time")
  
  # Plot detrended time series
  plot(detrended_ts, main = paste("Detrended", names(ts_list)[i]), 
       ylab = paste("Detrended", names(ts_list)[i]), xlab = "Time")
}
par(mfrow = c(1, 1))



par(mfrow = c(2, 2)) 
# Plot residuals for each time series
for (i in 1:length(residuals_list)) {
  ts_name <- names(residuals_list)[i]
  plot(residuals_list[[i]], type = "l", 
       main = paste("Residuals of", ts_name), 
       ylab = "Residuals", xlab = "Time")
}
par(mfrow = c(1, 1))


par(mfrow = c(2, 2))
for (i in seq_along(ts_list)) {
  original_ts <- ts_list[[i]]
  fit_trend <- fit_trend_list[[i]]
  
  # Create a time index for plotting (assuming monthly data)
  time_index <- 1:length(original_ts)
  
  # Plot original time series
  plot(time_index, original_ts, main = paste("Original", names(ts_list)[i]), 
       ylab = names(ts_list)[i], xlab = "Time", type = 'l', col = 'black', lwd = 2)
  
  # Add fitted trend
  lines(time_index, fit_trend, col = "red", lwd = 2)
}
par(mfrow = c(1, 1))


## ------------------------------------------ ##
# 2) Fit AR(1)
## ------------------------------------------ ##

ts_list <- detrended_ts_list

for (i in 1:length(ts_list)) {
  cat("Time Series", names(ts_list)[i], ":\n")
  
  # Fit the AR(1)
  AR1_model <- Arima(ts_list[[i]], order=c(1,0,0))
  # p = 1 = one autoregressive term
  # d = 0 = no differencing 
  # q = 0 = no moving average term
  
  #MLE
  #Arima (from forecast) and arima (from stats) functions are basically the same 
  
  AR1_sarima <- sarima(ts_list[[i]], 1, 0, 0) #same as Arima but summary gives more info
  
  # summary of fitted model
  print(summary(AR1_model))
  
  tsdiag(AR1_model) #diagnostic plots 
  
  # Extract residuals
  residuals <- residuals(AR1_model)
  
  # Plot residuals
  par(mfrow=c(2, 2))
  
  # Plot original time series and fitted values
  ts.plot(ts_list[[i]], main=paste("Original and Fitted for Time Series", 
                                   names(ts_list)[i]))
  lines(fitted(AR1_model), col="red")
  
  plot(residuals, main=paste("Residuals for Time Series", 
                             names(ts_list)[i]), ylab="Residuals")
  
  # Histogram of residuals
  hist(residuals, main=paste("Histogram of Residuals for Time Series", 
                             names(ts_list)[i]), xlab="Residuals")
  
  # Q-Q plot of residuals
  qqnorm(residuals)
  qqline(residuals, col="red")
  
  par(mfrow=c(2, 1)) 
  
  # Plot the ACF of the residuals
  #Acf(residuals, main=paste("ACF of Residuals for Time Series", 
  #                          names(ts_list)[i]))
  #acf1(residuals)
  acf2(residuals)
  
  #pacf(ts_list[[i]], main = paste("PACF of", ts_names[i]))
  
  # Perform the Ljung-Box test
  lb_test <- Box.test(residuals, lag=log(length(residuals)))
  print(lb_test)
  
  cat("\n")  
}
par(mfrow=c(1, 1))


par(mfrow=c(2, 1))
# Loop through each time series
for (i in 1:length(ts_list)) {
  
  # Plot the original time series
  ts.plot(ts_list[[i]], main = ts_names[i])
  
  # Fit the AR(1) model
  ar1_model <- Arima(ts_list[[i]], order=c(1,0,0))
  
  # Calculate fitted values
  ar1_fit <- fitted(ar1_model) #ar1_fit <- time_series[[i]] - residuals(ar1_model)
  
  # Add the fitted values to the plot
  points(ar1_fit, type = "l", col = 2, lty = 2)
  
  # Plot the residuals
  residuals_plot <- ts.plot(residuals(ar1_model), 
                            main=paste("Residuals for", 
                                       ts_names[i]))
}
par(mfrow=c(1, 1))


## ------------------------------------------ ##
# 3) Calculate AR coefficients
## ------------------------------------------ ##

# Create an empty data frame to store AR coefficients + other data
ar_info_df <- data.frame(ts_name = character(),
                         AR_coef = numeric(),
                         se = numeric(),
                         sigma2 = numeric(),
                         stringsAsFactors = FALSE)

# Loop through each time series
for (i in seq_along(ts_list)) {
  cat("Time Series", names(ts_list)[i], ":\n")
  
  # Fit the AR model with the best order
  best_model <- Arima(ts_list[[i]], order = c(1, 0, 0)) #MLE
  #Arima (from forecast) and arima (from stats) functions are basically the same 
  
  # Extract the AR coefficient
  ar_coef <- best_model$coef[1]
  se <- sqrt(diag(vcov(best_model)))[1]  # standard error
  sigma2 <- best_model$sigma2
  
  # Add information to the df
  ar_info_df <- rbind(ar_info_df, data.frame(ts_name = names(ts_list)[i],
                                             AR_coef = ar_coef,
                                             se = se,
                                             sigma2 = sigma2))
  
  cat("AR Coefficient:", ar_coef, "\n\n")
  cat("Standard Error (SE):", se, "\n")
  cat("Sigma squared (sigma2):", sigma2, "\n\n")
}

print(ar_info_df)

#write.csv(ar_info_df, "output/AR_coef_drivers_allAR1_afterDetrend.csv")

## ------------------------------------------ ##
###############################################################################
## ------------------------------------------ ##
###############################################################################

# THE END


## ------------------------------------------ ##
###############################################################################
## ------------------------------------------ ##
###############################################################################
## ------------------------------------------ ##
###############################################################################
## ------------------------------------------ ##
###############################################################################

### OLDER CODE
## ------------------------------------------ ##
# Check auto.arima 
## ------------------------------------------ ##

# list to store the ARIMA models
arima_models <- list()

# Loop through each time series and fit ARIMA model using auto arima
for (i in 1:length(ts_list)) {
  cat("Fitting ARIMA model to", names(ts_list)[i], ":\n")
  
  # Fit ARIMA model using auto.arima
  arima_models[[i]] <- auto.arima(ts_list[[i]], max.p = 10, 
                                  seasonal = FALSE, 
                                  ic = "aic") #here, got same results with aic and bic
  
  #summary of the fitted model
  print(summary(arima_models[[i]]))
  
  tsdiag(arima_models[[i]]) #diagnostic plots 
  
  residuals <- residuals(arima_models[[i]])
  
  # Plot residuals
  par(mfrow=c(2, 2))
  ts.plot(ts_list[[i]], main=paste("Original and Fitted for Time Series", 
                                   names(ts_list)[i]))
  lines(fitted(arima_models[[i]]), col="red")
  plot(residuals, main=paste("Residuals for Time Series", 
                             names(ts_list)[i]), ylab="Residuals")
  
  hist(residuals, main=paste("Histogram of Residuals for Time Series", 
                             names(ts_list)[i]), xlab="Residuals")
  qqnorm(residuals)
  qqline(residuals, col="red")
  
  par(mfrow=c(1, 1))
  Acf(residuals, main=paste("ACF of Residuals for Time Series", 
                            names(ts_list)[i]))
  
  # Ljung-Box test
  lb_test <- Box.test(residuals, lag=log(length(residuals)), type = c("Ljung-Box"))
  # I also looked at Box-Pierce results and very similar to Ljung-Box
  print(lb_test)
  
  cat("\n")
}



# df to store AR model information
ar_info_df2 <- data.frame(ts_name = character(),
                          AR_coef = numeric(),
                          ma1 = numeric(),
                          ma2 = numeric(),
                          se = numeric(),
                          sigma2 = numeric(),
                          stringsAsFactors = FALSE)

# Loop through each time series
for (i in seq_along(ts_list)) {
  cat("Time Series", names(ts_list)[i], ":\n")
  
  # Fit the AR model
  best_model <- auto.arima(ts_list[[i]], max.p = 10, seasonal = FALSE, ic = "aic")
  
  # Extract
  coef_names <- names(best_model$coef)
  ar_coef <- if ("ar1" %in% coef_names) best_model$coef["ar1"] else NA
  ma1 <- if ("ma1" %in% coef_names) best_model$coef["ma1"] else NA
  ma2 <- if ("ma2" %in% coef_names) best_model$coef["ma2"] else NA
  se <- sqrt(diag(vcov(best_model)))[1]  # Standard error
  sigma2 <- best_model$sigma2
  
  ar_info_df2 <- rbind(ar_info_df2, data.frame(ts_name = names(ts_list)[i],
                                               AR_coef = ar_coef,
                                               ma1 = ma1,
                                               ma2 = ma2,
                                               se = se,
                                               sigma2 = sigma2))
  
  cat("AR Coefficient:", ar_coef, "\n")
  cat("MA1 Coefficient:", ma1, "\n")
  cat("MA2 Coefficient:", ma2, "\n")
  cat("Standard Error (SE):", se, "\n")
  cat("Sigma squared (sigma2):", sigma2, "\n\n")
}


# Box.test(residuals(ar1_model), lag = 8, type = "Ljung-Box")
# doesnt matter what i do to the lag value; still significant
# small p value = suggest autocorrelation in residuals = model not great.?
library(readr)
library(purrr)

file_names <- list.files("output/", pattern = ".*Integrations.csv$", full.names = TRUE)

# Read all CSV files into a list of data frames
data_frames <- map(file_names, read_csv)

# unique data frame objects
for (i in seq_along(data_frames)) {
  # Extract unique identifier from the file name
  unique_id <- gsub(".*/(.*?)\\.csv", "\\1", file_names[i])
  
  var_name <- paste(unique_id, "data", sep = "_")
  
  # Assign data frame to unique variable name
  assign(var_name, data_frames[[i]])
}

rolling_variances <- zoo::rollapply(aoIntegrations_data$aoNorm, width = 12, FUN = var, align = "right", fill = NA)
# Plot the rolling variances over time
plot(rolling_variances, type = "l", 
     xlab = "Time", 
     ylab = "Variance", 
     main = "AO Rolling Variance")

# can add a horizontal line to visualize the average variance
abline(h = mean(rolling_variances, na.rm = TRUE), col = "red", lty = 2)


rolling_variances_amo <- zoo::rollapply(amoIntegrations_data$amoNorm, width = 12, FUN = var, align = "right", fill = NA)
# Plot the rolling variances over time
plot(rolling_variances_amo, type = "l", 
     xlab = "Time", 
     ylab = "Variance", 
     main = "AMO Rolling Variance")

# can add a horizontal line to visualize the average variance
abline(h = mean(rolling_variances_amo, na.rm = TRUE), col = "red", lty = 2)


rolling_variances_nao <- zoo::rollapply(naoIntegrations_data$naoNorm, width = 12, FUN = var, align = "right", fill = NA)
# Plot the rolling variances over time
plot(rolling_variances_nao, type = "l", 
     xlab = "Time", 
     ylab = "Variance", 
     main = "NAO Rolling Variance")

# can add a horizontal line to visualize the average variance
abline(h = mean(rolling_variances_nao, na.rm = TRUE), col = "red", lty = 2)


rolling_variances_gsi <- zoo::rollapply(gsiIntegrations_data$gsiNorm, width = 12, FUN = var, align = "right", fill = NA)
# Plot the rolling variances over time
plot(rolling_variances_gsi, type = "l", 
     xlab = "Time", 
     ylab = "Variance", 
     main = "GSI Rolling Variance")

# can add a horizontal line to visualize the average variance
abline(h = mean(rolling_variances_gsi, na.rm = TRUE), col = "red", lty = 2)


# variance 
#bartlett.test() or fligner.test()