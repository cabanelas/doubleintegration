################################################################################
#############          Pelagic Synthesis           #############################
#############             FEB-2025                 #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas 
################################################################################
##### CCE LTER
##### Sensitivity Analysis 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(tidyverse) #v2.0.0
library(doParallel) #better than parallel; v1.0.17
library(foreach)

## ------------------------------------------ ##
#            Biology Data -----
## ------------------------------------------ ##
bio <- read.csv(file.path("output",
                            "CCE",
                            "nsimplex_bio_output.csv")) %>%
  mutate(taxa = "Nsimplex") 

## ------------------------------------------ ##
#            PDO Data -----
## ------------------------------------------ ##
pdo <- read.csv(file.path("output",
                            "CCE",
                            "pdo_z_output.csv"))

## ------------------------------------------ ##
#            AR coefficients -----
## ------------------------------------------ ##
ar_bio <- read.csv(file.path("output",
                             "CCE",
                             "AR_coef_bio_CCE.csv"))

#ar_driver <- read.csv(file.path("output",
#                          "CCE",
#                          "AR_coef_PDOdriver_CCE.csv"))

# Extract AR(1) coefficients
#ar_coef_driver <- ar_driver$AR_coef[1]
ar_coef_bio <- ar_bio$AR_coef[1]

## ------------------------------------------ ##
#            Tidy Data -----
## ------------------------------------------ ##
sapply(bio, class)
sapply(pdo, class)

bio$date <- as.Date(bio$date)
pdo$time <- as.Date(pdo$time)

## ------------------------------------------ ##
#            Integrations -----
## ------------------------------------------ ##
# integration function (AR(1))
recursive_integration <- function(x, tau, dt = 1) {
  alpha <- 1 - dt / tau
  y <- numeric(length(x))
  y[1] <- x[1]
  for (i in 2:length(x)) {
    y[i] <- alpha * y[i - 1] + x[i] * dt
  }
  return(y)
}

# Z-score helpers
meanNaN <- function(x) mean(x, na.rm = TRUE)
sdNaN <- function(x) sd(x, na.rm = TRUE) #sample SD 

## ------------------------------------------ ##
#            Sensitivity -----
## ------------------------------------------ ##
taus <- seq(6, 120, by = 6)  # 0.5 to 10 years in 0.5-year steps (in months)

cor_vals <- numeric(length(taus))

for (i in seq_along(taus)) {
  tau <- taus[i]
  
  # Integrate PDO with current tau
  pdo_int <- recursive_integration(pdo$pdo_z, tau = tau, dt = 1)
  pdo_int_z <- (pdo_int - meanNaN(pdo_int)) / sdNaN(pdo_int)
  
  # Interpolate monthly PDO to yearly biology dates
  pdo_interp <- approx(pdo$time, pdo_int_z, xout = bio$date)$y
  
  # Correlate with annual biology anomalies
  cor_vals[i] <- cor(bio$Anomaly_yr, pdo_interp, use = "pairwise.complete.obs")
}

# Set up parallel cluster
num_cores <- parallel::detectCores() - 2
cl <- makeCluster(num_cores)
registerDoParallel(cl)


# Bootstrap settings
num_iterations <- 1000
n_driver <- nrow(pdo)
n_bio <- nrow(bio)
n_tau <- length(taus)

pdo_z <- pdo$pdo_z

# Run bootstrap loop
null_mat <- foreach(k = 1:num_iterations, .combine = rbind, .packages = "stats") %dopar% {
  # Simulate red-noise time series for driver and biology
  #red_driver <- arima.sim(n = n_driver, model = list(ar = ar_coef_driver))
  red_bio <- arima.sim(n = n_bio, model = list(ar = ar_coef_bio))
  
  result_row <- numeric(n_tau)
  
  for (j in seq_along(taus)) {
    tau <- taus[j]
    
    # Integrate red-noise driver
    driver_int <- recursive_integration(pdo_z, tau = tau, dt = 1)
    driver_int_z <- (driver_int - meanNaN(driver_int)) / sdNaN(driver_int)
    
    # Interpolate to match biology dates (assumes evenly spaced time points)
    #driver_interp <- approx(1:n_driver, driver_int_z, xout = seq_len(n_bio))$y
    driver_interp <- approx(pdo$time, driver_int_z, xout = bio$date)$y
    
    # Correlate integrated driver with red-noise biology
    result_row[j] <- cor(driver_interp, red_bio, use = "complete.obs")
  }
  
  result_row  # one row per iteration
}

# For each tau, get the 95th percentile of the null correlation distribution
sig_thresholds <- apply(null_mat, 2, quantile, probs = 0.95)

stopCluster(cl)

plot(taus / 12, cor_vals, type = "l", lwd = 2, col = "black",
     xlab = "Timescale τ (years)", ylab = "Correlation",
     main = "Sensitivity to τ (biological memory)",
     ylim = c(min(cor_vals, sig_thresholds) - 0.05, max(cor_vals, sig_thresholds) + 0.05))

# Add 95% null threshold (gray dashed line)
lines(taus / 12, sig_thresholds, lty = 2, col = "gray40")
