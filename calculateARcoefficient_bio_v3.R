################################################################################
#############          Pelagic Synthesis           #############################
#############             MAR-2024                 #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas
################################################################################
## Double Integration Analysis NES
# Script #3 : calculateARcoefficient_bio 

# script to detrend and calculate AR coefficients of biology time series

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
library(tibble) #v3.2.1
library(purrr) #v1.0.2
library(tidyr) #v1.3.1
library(astsa) #v2.1

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
bioTS <- read.csv(file.path("output",
                            "zscore_anomalies_3taxazp_28MAY.csv"))

bioTS <- bioTS %>%
  select(season, taxa, Region, Anomaly_yr, date)

#remove NAs
bioTS <- bioTS[complete.cases(bioTS$Anomaly_yr), ]

bioTS$date <- as.Date(bioTS$date, format="%Y-%m-%d")
bioTS$year <- format(bioTS$date, "%Y")
bioTS$year <- as.numeric(bioTS$year)

# Split - list of ts by season, taxa, and Region
bioTS_list <- split(bioTS, list(bioTS$season, bioTS$taxa, bioTS$Region))

# Remove empty elements
bioTS_list <- bioTS_list[sapply(bioTS_list, function(x) nrow(x) > 0)]

# names of the ts
# names == season.taxa.region
ts_names <- names(bioTS_list)

## ------------------------------------------ ##
#            Issue with missing years and making it a ts object -----
## ------------------------------------------ ##
unique_groups <- unique(bioTS[, c("season", "taxa", "Region")])

# Function to add missing years, fill with NA in Anomaly_yr column, and adjust date
add_missing_years <- function(group_df) {
  min_year <- min(group_df$year)
  max_year <- max(group_df$year)
  
  # Generate a sequence of all years from min to max
  all_years <- seq(min_year, max_year)
  
  # Create a data frame with all_years and all combinations of season, taxa, and Region
  # This ensures that all combinations of year, season, taxa, and Region are included
  all_combinations <- expand.grid(year = all_years,
                                  season = unique(group_df$season),
                                  taxa = unique(group_df$taxa),
                                  Region = unique(group_df$Region))
  
  # Merge with original data to ensure all years and combinations are included
  merged_data <- merge(all_combinations, group_df, 
                       by = c("year", "season", "taxa", "Region"), 
                       all.x = TRUE)
  
  # Sort by year if needed
  merged_data <- merged_data[order(merged_data$year), ]
  
  # Replace Anomaly_yr with NA where it's missing (assuming Anomaly_yr is numeric)
  merged_data$Anomaly_yr[is.na(merged_data$Anomaly_yr)] <- NA
  
  # Adjust date column by filling in missing years with a fabricated date
  # Assuming date is already of class Date, otherwise convert if needed
  if (!inherits(merged_data$date, "Date")) {
    merged_data$date <- as.Date(merged_data$date)
  }
  for (i in 2:nrow(merged_data)) {
    if (is.na(merged_data$Anomaly_yr[i])) {
      previous_date <- merged_data$date[i - 1]
      year <- merged_data$year[i]
      # Create a new date with the same month and day as the previous date
      new_date <- as.Date(paste(year, format(previous_date, "%m-%d"), sep = "-"))
      merged_data$date[i] <- new_date
    }
  }
  
  return(merged_data)
}

# Apply the function to each unique group combination
filled_data <- lapply(1:nrow(unique_groups), function(i) {
  group_df <- bioTS[bioTS$season == unique_groups$season[i] &
                      bioTS$taxa == unique_groups$taxa[i] &
                      bioTS$Region == unique_groups$Region[i], ]
  add_missing_years(group_df)
})

# Convert filled_data back to a single data frame if needed
filled_df <- do.call(rbind, filled_data)
bioTS_list1 <- split(filled_df, 
                     list(filled_df$season, 
                          filled_df$taxa, 
                          filled_df$Region))


## ------------------------------------------ ##
#            Make ts object -----
## ------------------------------------------ ##

convert_to_ts <- function(group_df) {
  # Subset relevant columns
  data_subset <- group_df[, c("year", "Anomaly_yr")]
  
  # Create ts object
  ts_data <- ts(data_subset$Anomaly_yr, start = min(data_subset$year), frequency = 1)
  
  return(ts_data)
}

# Apply the function to each group in bioTS_list
ts_list <- lapply(bioTS_list1, convert_to_ts)

bioTS_list <- ts_list

## ------------------------------------------ ##
# 1) Detrend time series ----> then AR(1)
## ------------------------------------------ ##

# Function to detrend time series using linear regression
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

detrended_results_list <- lapply(bioTS_list, detrend_time_series)

detrended_ts_list <- lapply(detrended_results_list, function(x) x$detrended_ts)
residuals_list <- lapply(detrended_results_list, function(x) x$residuals)
fit_trend_list <- lapply(detrended_results_list, function(x) x$fitted_trend)


detrended_df <- detrended_ts_list %>%
  map_dfr(~ as_tibble(as.numeric(.x)), .id = "series") %>%
  pivot_longer(cols = -series, names_to = "time", values_to = "value") %>%
  mutate(time = as.numeric(time))

# Function to extract date data from each time series object
extract_ts_data <- function(ts_obj) {
  tibble(
    year = time(ts_obj),
    value = as.numeric(ts_obj)
  )
}

# Apply the function to each element in the list and combine results
detrended_df <- detrended_ts_list %>%
  map_dfr(~ extract_ts_data(.x), .id = "series")

detrended_df <- detrended_df %>%
  separate(series, into = c("season", "taxa", "Region"), sep = "\\.")

# Merge detrended_df with bioTS based on season, taxa, Region, and year
merged_df <- detrended_df %>%
  left_join(bioTS %>% select(season, taxa, Region, year, date), 
            by = c("season", "taxa", "Region", "year"))

merged_df <- merged_df %>%
  rename(Anomaly_yr = value)

#write.csv(merged_df, "detrended_BIOLOGY_time_series.csv", row.names = TRUE)

## ------------------------------------------ ##
# PLOTS
## ------------------------------------------ ##

par(mfrow = c(2,1))
for (i in seq_along(bioTS_list)) {
  original_ts <- bioTS_list[[i]]
  detrended_ts <- detrended_ts_list[[i]]
  
  # Plot original time series
  plot(original_ts, main = paste("Original", names(bioTS_list)[i]), 
       ylab = names(bioTS_list)[i], xlab = "Time")
  
  # Plot detrended time series
  plot(detrended_ts, main = paste("Detrended", names(bioTS_list)[i]), 
       ylab = paste("Detrended", names(bioTS_list)[i]), xlab = "Time")
  
}
par(mfrow = c(1, 1))



par(mfrow = c(2, 2)) 
# Plot residuals for each time series
for (i in 1:length(residuals_list)) {
  ts_name <- names(residuals_list)[i] 
  plot(residuals_list[[i]], type = "l", main = paste("Residuals of", ts_name), ylab = "Residuals", xlab = "Time")
}
par(mfrow = c(1, 1))



par(mfrow = c(2, 3))
for (i in seq_along(bioTS_list)) {
  original_ts <- bioTS_list[[i]]
  fit_trend <- fit_trend_list[[i]]
  
  # Remove NA values from original_ts and corresponding fit_trend values
  na_mask <- !is.na(original_ts)
  original_ts <- original_ts[na_mask]
  fit_trend <- fit_trend[na_mask]
  
  # Create a time index for plotting
  time_index <- 1:length(original_ts)
  
  # Plot original time series
  plot(time_index, original_ts, main = paste(" ", names(bioTS_list)[i]), 
       ylab = names(bioTS_list)[i], xlab = "Time", type = 'l', col = 'blue', lwd = 2)
  
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
  
  AR1_sarima <- sarima(ts_list[[i]], 1, 0, 0)
  
  print(summary(AR1_model))
  
  tsdiag(AR1_model) #diagnostic plots 
  
  # Extract the residuals
  residuals <- residuals(AR1_model)
  
  # Plot the residuals
  par(mfrow=c(2, 2))  
  
  # Plot the original time series and the fitted values
  ts.plot(ts_list[[i]], main=paste("Original and Fitted for Time Series", 
                                   names(ts_list)[i]))
  lines(fitted(AR1_model), col="red")
  
  plot(residuals, main=paste("Residuals for Time Series", 
                             names(ts_list)[i]), ylab="Residuals")
  
  # Histogram of the residuals
  hist(residuals, main=paste("Histogram of Residuals for Time Series", 
                             names(ts_list)[i]), xlab="Residuals")
  
  # Q-Q plot of the residuals
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

#chart.ACFplus(ts)

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
  se <- sqrt(diag(vcov(best_model)))[1]  # Standard error
  sigma2 <- best_model$sigma2
  
  # Add the info to the df
  ar_info_df <- rbind(ar_info_df, data.frame(ts_name = names(ts_list)[i],
                                             AR_coef = ar_coef,
                                             se = se,
                                             sigma2 = sigma2))
  
  cat("AR Coefficient:", ar_coef, "\n\n")
  cat("Standard Error (SE):", se, "\n")
  cat("Sigma squared (sigma2):", sigma2, "\n\n")
}

print(ar_info_df)

# add sample size to each df (for building time series in bootstrap)
sample_sizes <- bioTS %>%
  group_by(season, taxa, Region) %>%
  summarise(sample_size = n()) %>%
  ungroup()

ar_info_df1 <- ar_info_df %>%
  mutate(season = sub("\\..*", "", ts_name),
         taxa = sub("^[^.]*\\.(.*?)\\..*", "\\1", ts_name),
         Region = sub("^[^.]*\\.[^.]*\\.(.*)", "\\1", ts_name))

ar_info_df1 <- ar_info_df1 %>%
  left_join(sample_sizes, by = c("season", "taxa", "Region"))

#write.csv(ar_info_df, "output/AR_coef_bio_allAR1_afterDetrend.csv")
#write.csv(ar_info_df1, "output/AR_coef_bio_allAR1_afterDetrend_wSampleSize.csv")

## ------------------------------------------ ##
###############################################################################
## ------------------------------------------ ##
###############################################################################

# THE END


