################################################################################
#############          Pelagic Synthesis           #############################
#############                 2024                 #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas 
################################################################################

## Double Integration Analysis NES
# Script #5 : Boot_parallel_custom_TAU

# script to generate two red-noise time series with the estimate AR coeff            

## STEP 2, 3, 4, & 5 of bootstrapping analysis
# 2. Generate surrogate red-noise time series with the same autoregression 
#coefficients
# 3. Calculate the correlation coefficient for each pair of time series
# 4. Repeat steps 2 and 3 for 10,000 iterations/Monte Carlo Simulation
# 5. Estimate the Probability Distribution Function (PDF) of corr coefficients

#showing the estimated PDF of the correlation coefficients
# Each plot based on the 10,000 correlation coefficients generated for that pair
# understand the distribution of the correlation coefficients and identify any 
#patterns or anomalies

######################################################################### 
# From Di Lorenzo and Ohman 

#The significance of the correlation coefficients is estimated from the 
#Probability Distribution Functions (PDFs) of the correlation coefficient 
#of two red-noise time series with the same autoregression coefficients as 
#estimated from the original signals. The PDFs are computed numerically by 
#generating 10,000 realizations of the correlation coefficient of two random 
#red-noise time series

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(dplyr) #v1.1.4
library(tidyr) #v1.3.1
library(here) #v1.0.1

#install.packages(c("progressr", "foreach", "doParallel"))
library(progress) #v1.2.3
library(foreach) #v1.5.2
library(doParallel) #better than parallel; v1.0.17

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##

AR_driver <- read.csv(file.path("output",
                                "AR_coef_drivers_allAR1_afterDetrend.csv"), 
                      header = T)

AR_biology <- read.csv(file.path("output",
                                 "AR_coef_bio_allAR1_afterDetrend_wSampleSize.csv"), 
                       header = T)

original_correlations <- read.csv(file.path("output",
                                            "IntCorrelations_DetrendBioDriver_6090240TAU.csv"), 
                                  header = T) %>%
  mutate(Pair_Name = paste(driver, season, taxa, region, sep = "."))


## ------------------------------------------ ##
#            singleIntCorrelations -----
## ------------------------------------------ ##
#
calculateIntegrations = function(data, tau, f = function(x) {mean(x, na.rm = T)}) {
  for (n in names(data)[-1]) {
    data[[paste0(n,'Norm')]] = (data[[n]] - mean(data[[n]], na.rm = T)) / sd(data[[n]]) #norm
    
    data[[paste0(n,'Int')]] = NA
    
    for (i in 2:nrow(data)) {
      k = data[,1] <= data[i,1] & data[,1] > data[i,1] - tau * 86400 #secs
      data[[paste0(n,'Int')]][i] = f(data[[paste0(n,'Norm')]][k])
    }
  }
  data
}

# Set tau values for different taxa
tau_values <- c(
  "ctyp" = 60, #2 months
  "calfin" = 240, #8 months
  "pseudo" = 90 #3 months
)

## ------------------------------------------ ##
#  PARAMETRIC BOOTSTRAPPING  ----
## ------------------------------------------ ##

#  Parallelize to speed up  ----

parallel::detectCores()
parallelly::availableCores()
numCores <- detectCores() - 2
cl <- makeCluster(numCores)
registerDoParallel(cl) 

# Preallocate lists
num_iterations <- 1000
num_drivers <- nrow(AR_driver)
num_biology <- nrow(AR_biology)
total_steps <- num_drivers * num_biology * num_iterations

# Initialize lists
original_cor_list <- vector("list", num_drivers * num_biology)
original_pval_list <- vector("list", num_drivers * num_biology)
cor_coeffs_list <- vector("list", num_drivers * num_biology)
ci_list <- vector("list", num_drivers * num_biology)
pvalue_list <- vector("list", num_drivers * num_biology)
pdfs <- vector("list", num_drivers * num_biology)
results <- vector("list", num_drivers * num_biology)

# Function to generate red-noise time series
generateRedNoise <- function(n, ar_coefs) {
  arima.sim(n = n, model = list(ar = ar_coefs))
}

# Function to interpolate time series
interpolateDriver <- function(driver_time, driver_data, abundance_time) {
  interpolated_driver <- approx(driver_time, driver_data, xout = abundance_time)$y
  return(interpolated_driver)
}

# Main loop with parallel processing
counter <- 1
start_time <- Sys.time()

set.seed(123)

results <- foreach(i = 1:num_drivers, 
                   .combine = 'rbind', 
                   .packages = c("progress", "foreach", "doParallel")) %dopar% {
  ar_coefs_driver <- na.omit(as.numeric(AR_driver[i, 3]))
  driver_results <- list()
  
  for (j in 1:num_biology) {
    ar_coefs_bio <- na.omit(as.numeric(AR_biology[j, 3]))
    sample_size_bio <- AR_biology$sample_size[j]
    
    ts_name <- AR_biology$ts_name[j]
    taxa <- strsplit(ts_name, "\\.")[[1]][2]
    tau <- tau_values[taxa]
    
    cor_coeffs_pair <- numeric(num_iterations)
    cor_coeffs_pair_integrated <- numeric(num_iterations)
    
    for (k in 1:num_iterations) {
      red_noise_driver <- generateRedNoise(610, ar_coefs_driver)
      red_noise_bio <- generateRedNoise(sample_size_bio, ar_coefs_bio)
      
      data_driver <- data.frame(time = 1:610, driver = red_noise_driver)
      
      integrated_drivers <- lapply(tau_values, 
                                   function(tau) calculateIntegrations(data_driver, tau))
      
      integrated_driver <- integrated_drivers[[taxa]]$driverInt
      interpolated_driver <- interpolateDriver(1:610, 
                                               red_noise_driver, 
                                               1:sample_size_bio)
      interpolated_integrated_driver <- interpolateDriver(1:610, 
                                                          integrated_driver, 
                                                          1:sample_size_bio)
      
      cor_coeffs_pair[k] <- cor(interpolated_driver, 
                                red_noise_bio, 
                                use = "complete.obs")
      cor_coeffs_pair_integrated[k] <- cor(interpolated_integrated_driver, 
                                           red_noise_bio, use = "complete.obs", 
                                           method = "pearson")
      
      if (k %% 200 == 0) {
        write(paste("Driver:", i, "Group:", j, "Iteration:", k), "counter.txt")
      }
    }
    
    driver_name <- tolower(sub("ts", "", AR_driver$ts_name[i]))
    group_name <- sub("ts", "", AR_biology$ts_name[j])
    driver_parts <- unlist(strsplit(driver_name, "\\."))
    group_parts <- unlist(strsplit(group_name, "\\."))
    
    standardized_pair_name <- paste(driver_parts[1], 
                                    group_parts[3], 
                                    group_parts[1], 
                                    group_parts[2], 
                                    sep = ".")
    standardized_pair_name_integrated <- paste(standardized_pair_name, 
                                               "integrated", sep = ".")
    
    original_cor_norm <- original_correlations$cor_Norm[which(original_correlations$Pair_Name == standardized_pair_name)]
    original_cor_int <- original_correlations$cor_Int[which(original_correlations$Pair_Name == standardized_pair_name)]
    
    original_pval_norm <- original_correlations$pval_Norm[which(original_correlations$Pair_Name == standardized_pair_name)]
    original_pval_int <- original_correlations$pval_Int[which(original_correlations$Pair_Name == standardized_pair_name)]
    
    cor_coeffs_pair_mean <- mean(cor_coeffs_pair)
    cor_coeffs_pair_integrated_mean <- mean(cor_coeffs_pair_integrated)
    
    ci_pair <- quantile(cor_coeffs_pair, 
                        probs = c(0.025, 0.975))
    ci_pair_integrated <- quantile(cor_coeffs_pair_integrated, 
                                   probs = c(0.025, 0.975))
    
    pval_pair <- mean(abs(cor_coeffs_pair) >= abs(original_cor_norm))
    pval_pair_integrated <- mean(abs(cor_coeffs_pair_integrated) >= abs(original_cor_int))
    
    driver_results[[j]] <- data.frame(
      driver_name = driver_name,
      group_name = group_name,
      standardized_pair_name = standardized_pair_name,
      original_cor_norm = original_cor_norm,
      original_cor_int = original_cor_int,
      original_pval_norm = original_pval_norm,
      original_pval_int = original_pval_int,
      cor_coeffs_pair_mean = cor_coeffs_pair_mean,
      cor_coeffs_pair_integrated_mean = cor_coeffs_pair_integrated_mean,
      ci_pair_lower = ci_pair[1],
      ci_pair_upper = ci_pair[2],
      ci_pair_integrated_lower = ci_pair_integrated[1],
      ci_pair_integrated_upper = ci_pair_integrated[2],
      pval_pair = pval_pair,
      pval_pair_integrated = pval_pair_integrated
    )

  }
  
  do.call(rbind, driver_results)
}

# Stop the parallel cluster
stopCluster(cl)

end_time <- Sys.time()
end_time - start_time
cat("Total computation time: ", end_time - start_time, "\n")
#Total computation time:  5.151053 (initial run, before streamlining code)
#Total computation time:  1.161409 (with just 2 taxa)
#Total computation time:  2.546351 (a bit less with all 3 taxa)

results_df <- results
check_significance <- function(results_df) {
  results_df$significant_original <- with(results_df, 
                                          original_cor_norm < ci_pair_lower | 
                                            original_cor_norm > ci_pair_upper)
  
  results_df$significant_integrated <- with(results_df, 
                                            original_cor_int < ci_pair_integrated_lower | 
                                            original_cor_int > ci_pair_integrated_upper)
  return(results_df)
}

# Apply the function
significance_results <- check_significance(results_df)



results_df5 <- results_df %>%
  mutate(
    original_cor_norm = as.numeric(original_cor_norm),
    original_cor_int = as.numeric(original_cor_int),
    pval_pair = as.numeric(pval_pair),
    pval_pair_integrated = as.numeric(pval_pair_integrated),
    significance_percentage_norm = (1 - pval_pair) * 100,
    significance_percentage_integrated = (1 - pval_pair_integrated) * 100
  )


results_4 <- significance_results %>%
  mutate(
    original_cor_norm = as.numeric(original_cor_norm),
    original_cor_int = as.numeric(original_cor_int),
    pval_pair = as.numeric(pval_pair),
    pval_pair_integrated = as.numeric(pval_pair_integrated),
    significance_percentage_norm = (1 - pval_pair) * 100,
    significance_percentage_integrated = (1 - pval_pair_integrated) * 100
  )
#write.csv(results_4, "results_boot_6090240TAU_RAW.csv")
