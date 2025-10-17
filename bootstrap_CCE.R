################################################################################
#############          Pelagic Synthesis           #############################
#############                 2024                 #############################
#############         CCE - Double Integration     #############################
## by: Alexandra Cabanelas 
################################################################################
## CCE LTER
## Driver = Pacific Decadal Oscillation (PDO)
## Biology = Nyctiphanes simplex (abundance m2)
## season = spring (Feb, Mar, Apr, May)
## years = 1951 - 2021
## night tows only

# Script #4 : boostrap_CCE

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
#library(dplyr) #v1.1.4
#library(tidyr) #v1.3.1
#library(here) #v1.0.1

#install.packages(c("foreach", "doParallel"))
library(foreach) #v1.5.2
library(doParallel) #better than parallel; v1.0.17
library(tidyverse)
library(cowplot) #v1.1.3
library(doRNG) # so i can parallel while also reproducible

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
AR_driver <- read.csv(file.path("output",
                                "CCE",
                                "ARcoef_PDOdriver_CCE.csv"), 
                      header = T) # created in script#1

AR_biology <- read.csv(file.path("output",
                                 "CCE",
                                 "ARcoef_Nsimplexbio_CCE.csv"), 
                       header = T) #created in script#2

original_correlations <- read.csv(file.path("output",
                                            "CCE",
                                            "cor_Nsimplex_CCE.csv"), 
                                  header = T) #created in script#3
# for plots..
PDO <- read.csv(file.path("output",
                          "CCE",
                          "PDOint_Nsimplex_CCE.csv"), 
                header = T) %>%
  mutate(time = as.Date(time, 
                        format = "%Y-%m-%d")) #created in script#3

euphs <- read.csv(file.path("output",
                            "CCE",
                            "nsimplex_anomalies.csv")) %>% #created in script#3
  mutate(date = as.Date(date))

## ------------------------------------------ ##
#     Function for integrations -----
## ------------------------------------------ ##
recursive_integration <- function(x, tau, dt = 1) {
  # calculate autoregressive weight (how much past value is retained)
  # tau = memory in months, dt = time step 1 month
  alpha <- 1 - dt / tau # autoregressive weight; alpha bio, Euler
  y <- numeric(length(x)) #bio response
  y[1] <- x[1] # initial condition; (MATLAB: o(1) = forc.sig(1))
  for (i in 1:(length(x) - 1)) {
    y[i + 1] <- alpha * y[i] + x[i] * dt #forward recursion loop
  }
  return(y)
}

## ------------------------------------------ ##
#  PARAMETRIC BOOTSTRAPPING  ----
## ------------------------------------------ ##

#  --- Parallelize to speed up -----
parallel::detectCores()
parallelly::availableCores()
numCores <- detectCores() - 2
cl <- makeCluster(numCores)
registerDoParallel(cl) 

set.seed(123)  # ensure reproducibility

#  --- define AR(1) coefficients and sample sizes -----
ar_coef_driver <- as.numeric(AR_driver$AR_coef)  #same as AR_driver$AR_coef[1] 
ar_coef_bio <- as.numeric(AR_biology$AR_coef)  
n_driver <- AR_driver$n #993  #length of PDO starting at 1941 **CHECK
n_bio <- AR_biology$n  # length of biological time series (includig NAs)

#  --- define parameters -----
num_iterations <- 10000
tau_bio <- 24  # biological memory in months
dt <- 1        # monthly timestep

## ------------------------------------------ ##
#  Bootstrap Loop  ----
## ------------------------------------------ ##
start_time <- Sys.time()

cor_results <- foreach(k = 1:num_iterations, .combine = rbind, .packages = "stats") %dorng% {
  
  # generate new red-noise time series each iteration
  red_noise_driver_k <- arima.sim(n = n_driver, 
                                  model = list(ar = ar_coef_driver))
  red_noise_bio_k    <- arima.sim(n = n_bio, 
                                  model = list(ar = ar_coef_bio))
  
  # Normalize driver 
  driver_norm <- scale(red_noise_driver_k)[, 1] #or as.numeric(scale(red_noise_driver_k))
  
  # Apply integration
  driver_int <- recursive_integration(driver_norm, tau = tau_bio, dt = dt)
  # Z-score after integration
  driver_int_z <- scale(driver_int)[, 1]
  
  # interpolate driver to match bio time points 
  idx_yearly <- seq(1, n_driver, length.out = n_bio)
  interpolated_driver <- approx(1:n_driver, driver_norm, xout = idx_yearly)$y
  interpolated_int    <- approx(1:n_driver, driver_int_z, xout = idx_yearly)$y

  # Compute correlations
  data.frame(
    cor_Norm = cor(interpolated_driver, red_noise_bio_k, use = "complete.obs"),
    cor_Int = cor(interpolated_int, red_noise_bio_k, use = "complete.obs")
  )
}

# stop the parallel cluster
stopCluster(cl)

end_time <- Sys.time()
end_time - start_time
cat("Total computation time: ", end_time - start_time, "\n")
#Total computation time: 2.94 mins; 2.66 mins; 6.86 secs

## ------------------------------------------ ##
#  Summarize Bootstrap Results
## ------------------------------------------ ##
final_results <- data.frame(
  # original correlations & (conventional) pvalues
  original_cor_norm = original_correlations$cor_Norm,
  original_cor_int = original_correlations$cor_Int,
  original_pval_norm = original_correlations$pval_Norm,
  original_pval_int = original_correlations$pval_Int,
  
  # mean bootstrapped correlations; near 0 = expected under null
  cor_coeffs_pair_mean = mean(cor_results$cor_Norm),
  cor_coeffs_pair_integrated_mean = mean(cor_results$cor_Int),
  
  # 95% confidence intervals of null distribution
  ci_pair_lower = quantile(cor_results$cor_Norm, probs = 0.025),
  ci_pair_upper = quantile(cor_results$cor_Norm, probs = 0.975),
  ci_pair_integrated_lower = quantile(cor_results$cor_Int, probs = 0.025),
  ci_pair_integrated_upper = quantile(cor_results$cor_Int, probs = 0.975),
  
  # Empirical/Monte Carlo p-values
  pval_pair = mean(abs(cor_results$cor_Norm) >= abs(original_correlations$cor_Norm)),
  pval_pair_integrated = mean(abs(cor_results$cor_Int) >= abs(original_correlations$cor_Int))
)

# Flag significance based on bootstrap confidence intervals
final_results <- final_results %>%
  mutate(
    significant_original = original_cor_norm < ci_pair_lower | 
      original_cor_norm > ci_pair_upper,
    significant_integrated = original_cor_int < ci_pair_integrated_lower | 
      original_cor_int > ci_pair_integrated_upper
  )
final_results
#write.csv(final_results, "output/CCE/boot_results_NsimplexCCE.csv")

## ------------------------------------------ ##
# Sensitivity Analysis: Vary τ and compute correlation
## ------------------------------------------ ##

# range of tau values (in months)
tau_range <- seq(6, 120, by = 3)  # 0.5-10 years in 3month steps

#initialize results df
tau_sensitivity <- data.frame(tau = tau_range, cor = NA) 
#cor_vals <- numeric(length(tau_range))

# Loop over tau values
for (i in seq_along(tau_range)) {
  
  tau_i <- tau_range[i]
  
  # integrate PDO
  pdo_int_i <- recursive_integration(PDO$pdo_z, tau = tau_i, dt = 1)
  # z-score
  pdo_int_z_i <- as.numeric(scale(pdo_int_i)) 
  
  # Interpolate to biological sampling dates
  PDO_raw_interp <- approx(PDO$time, PDO$pdo_z, xout = euphs$date, rule = 1)$y
  PDO_int_interp <- approx(PDO$time, pdo_int_z_i, xout = euphs$date, rule = 1)$y
  
  # Compute correlation with biological anomalies
  cor_raw <- cor(PDO_raw_interp, euphs$Anomaly_yr, use = "complete.obs")
  cor_i <- cor(PDO_int_interp, euphs$Anomaly_yr, use = "complete.obs")
  
  # Store result
  tau_sensitivity$cor[i] <- cor_i
}

plot(tau_range / 12, tau_sensitivity$cor, type = "l", lwd = 2, col = "black", 
     xlab = "Timescale τ (years)", 
     ylab = "Correlation",
     main = "Sensitivity to τ (biological memory)")

library(signal)
lowpass_sensitivity <- data.frame(tau = tau_range, cor = NA)

for (i in seq_along(tau_range)) {
  tau_i <- tau_range[i]
  cutoff_freq <- 1 / tau_i  # cycles per month
  
  # Butterworth low-pass filter
  bf <- butter(n = 2, W = cutoff_freq, type = "low", plane = "z")
  pdo_filtered <- filtfilt(bf, PDO$pdo_z)
  pdo_filtered_z <- as.numeric(scale(pdo_filtered))
  
  PDO_interp <- approx(PDO$time, pdo_filtered_z, xout = euphs$date, rule = 2)$y
  cor_i <- cor(PDO_interp, euphs$Anomaly_yr, use = "complete.obs")
  lowpass_sensitivity$cor[i] <- cor_i
}

# convert tau to years
lowpass_sensitivity$tau_years <- tau_range / 12
tau_sensitivity$tau_years <- tau_sensitivity$tau / 12

#  results
print(tau_sensitivity)
cor_threshold_95 <- quantile(cor_results$cor_Int, probs = 0.95)

combined_df <- rbind(
  data.frame(tau = tau_sensitivity$tau_years, cor = tau_sensitivity$cor, Type = "Integrated PDO"),
  data.frame(tau = lowpass_sensitivity$tau_years, cor = lowpass_sensitivity$cor, Type = "Low-pass PDO")
)

ggplot(tau_sensitivity, aes(x = tau_years, y = cor)) +
  geom_hline(yintercept = cor_threshold_95, linetype = "dashed", 
             color = "gray40") +
  annotate("rect", xmin = 0.5, xmax = 2, ymin = -Inf, ymax = Inf, 
           alpha = 0.5, fill = "gray") +
  geom_line(color = "black", linewidth = 1.2) +
  geom_point(color = "blue", size = 2) +
  geom_hline(yintercept = cor_raw, color = "red", linewidth = 1.2) +
  labs(title = expression(italic(N.*simplex)~"sensitivity to biological damping timescale"),
       subtitle = paste("Black = integrated PDO; Red = raw PDO\nDashed line = 95% significance (r =", round(cor_threshold_95, 2), ")"),
       x = expression(paste("Biological memory ", tau[bio], " [years]")),
       y = "Correlation coefficient (r)",
       caption = "CCE – spring") +
  theme_minimal(base_size = 14)

peak <- tau_sensitivity[which.max(tau_sensitivity$cor), ]
ggplot(combined_df, aes(x = tau, y = cor, color = Type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = cor_threshold_95, 
             linetype = "dashed", color = "gray40") +
  annotate("rect", xmin = 0.5, xmax = 2, ymin = -Inf, ymax = Inf, 
           alpha = 0.3, fill = "gray") +
  geom_point(data = peak, aes(x = tau_years, y = cor), shape = 8, size = 3, color = "black") + 
  scale_color_manual(values = c("Integrated PDO" = "black", 
                                "Low-pass PDO" = "red")) +
  labs(title = expression(italic(N.*simplex)~"sensitivity to biological damping timescale"),
       subtitle = paste("Dashed line = 95% significance (r =", round(cor_threshold_95, 2), ")"),
       x = expression(paste("Timescale ", tau, " [years]")),
       y = "Correlation coefficient (r)",
       caption = "CCE – spring") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  scale_x_continuous(
    breaks = seq(0, 10, by = 2),
    limits = c(0, 10)
  ) 

## ------------------------------------------ ##
#     Plot PDFs -----
## ------------------------------------------ ##
obs_norm <- original_correlations$cor_Norm
obs_int <- original_correlations$cor_Int

cor_df <- cor_results %>%
  rename(Original = cor_Norm,
         Integrated = cor_Int) %>%
  pivot_longer(cols = c(Original, Integrated), 
               names_to = "Type", values_to = "Correlation") %>%
  mutate(Type = factor(Type, levels = c("Original", "Integrated")))

obs_df <- data.frame(
  Type = factor(c("Original", "Integrated"), 
                levels = c("Original", "Integrated")),
  Obs = c(obs_norm, obs_int)
)

(PDF_plot <- ggplot(cor_df, aes(x = Correlation, fill = Type)) +
    geom_density(alpha = 0.5, color = NA) +
    geom_vline(data = obs_df, aes(xintercept = Obs),
               color = "black", linetype = "dashed", linewidth = 1) +
    facet_wrap(~Type, scales = "free", ncol = 2) + 
    theme_minimal(base_size = 14) +
    theme(legend.position = "none") +
    labs(title = "Bootstrapped Correlation Distributions", #PDF
         subtitle = "Dashed line = observed correlation",
         x = "Correlation coefficient (r)",
         y = "Probability density",
         caption = "CCE – N. simplex – spring")
)

# export and save plot
#ggsave("figures/CCE/bootPDF_CCE_Ns_spring.png", PDF_plot, 
#       width = 10, height = 6, dpi = 300, bg = "white")

#shape and spread of both null distributions
plot(density(cor_results$cor_Norm), col = "red", lwd = 2,
     main = "Density of Bootstrapped Correlations",
     xlab = "Correlation (r)", ylim = c(0, max(density(cor_results$cor_Norm)$y)))
lines(density(cor_results$cor_Int), col = "blue", lwd = 2)
legend("topright", legend = c("PDO", "Integrated PDO"), col = c("red", "blue"), lwd = 2)
ci <- quantile(cor_results$cor_Int, probs = c(0.025, 0.975))
abline(v = ci, col = "darkgreen", lty = 2)

cor_long <- cor_results %>%
  pivot_longer(cols = c(cor_Norm, cor_Int),
               names_to = "Type", values_to = "r")

ggplot(cor_long, aes(x = r, fill = Type)) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = final_results$original_cor_norm),
             color = "blue", linetype = "solid") +
  geom_vline(aes(xintercept = final_results$original_cor_int),
             color = "blue", linetype = "dotted") +
  labs(
    title = "Comparison of null distributions",
    x = "Correlation coefficient (r)",
    y = "Density"
  ) +
  scale_fill_manual(values = c("cor_Norm" = "#999999", "cor_Int" = "#1B9E77"),
                    labels = c("Raw", "Integrated")) +
  theme_minimal(base_size = 14)

boxplot(cor_results$cor_Norm, cor_results$cor_Int,
        names = c("PDO", "Integrated PDO"),
        col = c("lightblue", "lightgreen"),
        main = "Bootstrap Correlation Distributions",
        ylab = "Correlation (r)")

## ------------------------------------------ ##
#     Plots -----
## ------------------------------------------ ##'

cor_CCE <- original_correlations

cor_norm <- cor_CCE$cor_Norm[1]
cor_int <- cor_CCE$cor_Int[1]

pval_norm <- final_results$pval_pair[1]
pval_int <- final_results$pval_pair_integrated[1]

(NormPlot <- ggplot() +
    # Driver time series (PDO normalized)
    geom_line(data = PDO, aes(x = time, y = pdo_z, color = "PDO"), 
              linewidth = 1.2) +
    # Nsimplex time series (Abundance)
    geom_line(data = euphs, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
              linewidth = 1.2) +
    # annotate correlation coefficient
    annotate("text", x = as.Date("1960-01-01"), y = 3,
             label = sprintf("r == %.2f", cor_norm),
             parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
    # annotate bootstrapped p-value
    annotate("text", x = as.Date("1960-01-01"), y = 2.4,
             label = if (pval_norm < 0.001) {
               as.expression(bquote(italic(p) < 0.001))
             } else {
               as.expression(bquote(italic(p) == .(round(pval_norm, 3))))
             },
             parse = TRUE, hjust = 0, vjust = 1, size = 5.5) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 19, color = "black", face = "bold"),
          axis.text.y = element_text(size = 15, color = "black"),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 15, face = "bold"),
          legend.key.size = unit(2, "lines"),
          legend.position = "bottom",
          plot.title = element_blank()) +
    # x-axis settings
    scale_x_date(breaks = seq(as.Date("1950-01-01"), as.Date("2022-01-01"), 
                              by = "10 years"),
                 date_labels = "%Y", expand = c(0, 0)) +
    coord_cartesian(xlim = as.Date(c("1950-01-01", "2022-01-01"))) +
    # manual color legend
    scale_color_manual(values = c("PDO" = "red", "Abundance" = "blue")) +
    labs(y = "PDO")
)

(Int1Plot <- ggplot() +
    # Driver time series (Integrated PDO)
    geom_line(data = PDO, aes(x = time, y = pdo_int_z, color = "PDO"), 
              linewidth = 1.2) +
    # Nsimplex time series (Abundance)
    geom_line(data = euphs, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
              linewidth = 1.2) +
    # annotate correlation coefficient
    annotate("text", x = as.Date("1960-01-01"), y = 3,
             label = sprintf("r == %.2f", cor_int),
             parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
    # annotate bootstrapped p-value
    annotate("text", x = as.Date("1960-01-01"), y = 2.4,
             label = if (pval_int < 0.001) {
               as.expression(bquote(italic(p) < 0.001))
             } else {
               as.expression(bquote(italic(p) == .(round(pval_int, 3))))
             },
             parse = TRUE, hjust = 0, vjust = 1, size = 5.5) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 19, color = "black", face = "bold"),
          axis.text = element_text(size = 15, color = "black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 15, face = "bold"),
          legend.key.size = unit(2, "lines"),
          plot.title = element_blank()) +
    # x-axis settings
    scale_x_date(breaks = seq(as.Date("1950-01-01"), as.Date("2022-01-01"), 
                              by = "10 years"),
                 date_labels = "%Y", expand = c(0, 0)) +
    coord_cartesian(xlim = as.Date(c("1950-01-01", "2022-01-01"))) +
    # manual color legend
    scale_color_manual(values = c("PDO" = "red", "Abundance" = "blue")) +
    # y-axis label
    labs(y = "Integrated PDO")
)

legend <- cowplot::get_plot_component(NormPlot, 'guide-box-bottom', 
                                      return_all = TRUE)
cowplot::ggdraw(legend)

# remove legend
NormPlot <- NormPlot + theme(legend.position = "none")
Int1Plot <- Int1Plot + theme(legend.position = "none")

# combine plots into one grid
combined_plots <- plot_grid(NormPlot, Int1Plot,
                            ncol = 1, align = "v")

# add plot title 
title <- ggdraw() +
  draw_label(label = bquote("CCE " * italic("Nyctiphanes simplex") * " Spring"), 
             size = 18) +
  theme(plot.background = element_rect(fill = "transparent"))


# combine title and plots
(final_plot <- plot_grid(title, combined_plots, legend, ncol = 1, 
                         rel_heights = c(0.07, 0.8, 0.07))) 

# export and save plot
#ggsave("figures/CCE/boot_CCE_Ns_spring.png", 
#       final_plot, 
#       width = 8, height = 8, dpi = 300, 
#       bg = "white")




####
run_bootstrap <- function(ar_coef_driver, ar_coef_bio, n_driver, n_bio, tau_bio, num_iterations = 10000) {
  # bootstrap logic here
  return(final_results)
}
#####

## extra PDF plot
summary_df <- data.frame(
  Type = c("Raw", "Integrated"),
  Observed = c(final_results$original_cor_norm, final_results$original_cor_int),
  MeanNull = c(final_results$cor_coeffs_pair_mean, final_results$cor_coeffs_pair_integrated_mean),
  CI_lower = c(final_results$ci_pair_lower, final_results$ci_pair_integrated_lower),
  CI_upper = c(final_results$ci_pair_upper, final_results$ci_pair_integrated_upper)
)

ggplot(summary_df, aes(x = Type)) +
  geom_pointrange(aes(y = MeanNull, ymin = CI_lower, ymax = CI_upper),
                  color = "grey50", size = 1.2) +
  geom_point(aes(y = Observed), color = "blue", size = 3) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    title = "Observed vs. null mean correlations",
    y = "Correlation coefficient (r)"
  ) +
  theme_minimal(base_size = 14)


