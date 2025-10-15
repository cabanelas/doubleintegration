################################################################################
#############          Pelagic Synthesis           #############################
#############                 2024                 #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas 
################################################################################

## Double Integration Analysis NGA
# Script #5 : boostrap_NGA

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

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
AR_driver <- read.csv(file.path("output",
                                "NGA",
                                "AR_coef_PDOdriver_NGA.csv"), 
                      header = T)

AR_biology <- read.csv(file.path("output",
                                 "NGA",
                                 "AR_coef_NeoCriSpbio_NGA.csv"), 
                       header = T) %>%
  mutate(n = 25) 

original_correlations <- read.csv(file.path("output",
                                            "NGA",
                                            "weightedCor_NGA.csv"), 
                                  header = T) #CHECK VERSION OF THIS FILE TO USE

## ------------------------------------------ ##
#     Function for weighted integrations -----
## ------------------------------------------ ##
calculateIntegrations = function(data, tau = 180, 
                                 f = function(x, w) {sum(x * w, na.rm = T) / sum(w, na.rm = T)}) {
  for (n in names(data)[-1]) {
    data[[paste0(n,'Norm')]] = (data[[n]] - mean(data[[n]], na.rm = T)) / sd(data[[n]]) # norm
    
    data[[paste0(n,'Int')]] = NA
    #data[[paste0(n,'DInt')]] = NA
    
    for (i in 2:nrow(data)) {
      #first integration
      k1 = data[,1] <= data[i,1] & data[,1] > data[i,1] - tau * 86400 # secs
      tdiff1 = as.numeric(data[i,1] - data[k1,1]) # convert difftime to numeric
      w1 = exp(-tdiff1 / (tau * 86400 / log(2))) # weights
      data[[paste0(n,'Int')]][i] = f(data[[paste0(n,'Norm')]][k1], w1)
      
      # second integration
      #k2 = data[,1] <= data[i,1] & data[,1] > data[i,1] - tau * 86400
      #tdiff2 = as.numeric(data[i,1] - data[k2,1])
      #w2 = exp(-tdiff2 / (tau * 86400 / log(2)))
      #data[[paste0(n, 'DInt')]][i] = f(data[[paste0(n, 'Int')]][k2], w2)
    }
 }
  data
}

## ------------------------------------------ ##
#  PARAMETRIC BOOTSTRAPPING  ----
## ------------------------------------------ ##

#  Parallelize to speed up  ----

parallel::detectCores()
parallelly::availableCores()
numCores <- detectCores() - 2
cl <- makeCluster(numCores)
registerDoParallel(cl) 
start_time <- Sys.time()

set.seed(123)  # ensure reproducibility

# define AR(1) coefficient and sample size
ar_coef_driver <- as.numeric(AR_driver$AR_coef)  
ar_coef_bio <- as.numeric(AR_biology$AR_coef)  
n_driver <- 993  #length of PDO starting at 1941
n_bio <- AR_biology$n  # length of biological time series 
num_iterations <- 10000

cor_results <- foreach(k = 1:num_iterations, .combine = rbind, .packages = "stats") %dopar% {
  
  # generate new red-noise time series each iteration
  red_noise_driver_k <- arima.sim(n = n_driver, model = list(ar = ar_coef_driver))
  red_noise_bio_k <- arima.sim(n = n_bio, model = list(ar = ar_coef_bio))
  
  # create driver data frame (just indices)
  driver_data <- data.frame(time = 1:n_driver, driver = red_noise_driver_k)
  
  # calculate integrations 
  driverInt <- calculateIntegrations(driver_data, tau = 365 * 1)
  
  # extract relevant columns
  driver1 <- driverInt[, c("time", "driverNorm")]
  int1 <- driverInt[, c("time", "driverInt")]
  
  # interpolate driver to match bio time points 
  interpolated_driver <- approx(driver1$time, driver1$driverNorm, xout = 1:n_bio)$y
  interpolated_int1 <- approx(int1$time, int1$driverInt, xout = 1:n_bio)$y
  
  # compute correlations
  data.frame(
    cor_Norm = cor(interpolated_driver, red_noise_bio_k, use = "complete.obs"),
    cor_Int = cor(interpolated_int1, red_noise_bio_k, use = "complete.obs")
  )
}
# stop the parallel cluster
stopCluster(cl)

end_time <- Sys.time()
end_time - start_time
cat("Total computation time: ", end_time - start_time, "\n")
#Total computation time: 2.71 mins

final_results <- data.frame(
  # original correlations 
  original_cor_norm = original_correlations$cor_Norm,
  original_cor_int = original_correlations$cor_Int,
  # original (conventional) pvalues
  original_pval_norm = original_correlations$pval_Norm,
  original_pval_int = original_correlations$pval_Int,
  # mean bootstrapped correlations; near 0 = expected under null
  cor_coeffs_pair_mean = mean(cor_results$cor_Norm),
  cor_coeffs_pair_integrated_mean = mean(cor_results$cor_Int),
  # 95% confidence intervals of distribution
  ci_pair_lower = quantile(cor_results$cor_Norm, probs = 0.025),
  ci_pair_upper = quantile(cor_results$cor_Norm, probs = 0.975),
  ci_pair_integrated_lower = quantile(cor_results$cor_Int, probs = 0.025),
  ci_pair_integrated_upper = quantile(cor_results$cor_Int, probs = 0.975),
  pval_pair = mean(abs(cor_results$cor_Norm) >= abs(original_correlations$cor_Norm)),
  pval_pair_integrated = mean(abs(cor_results$cor_Int) >= abs(original_correlations$cor_Int))
)

final_results <- final_results %>%
  mutate(
    significant_original = original_cor_norm < ci_pair_lower | 
      original_cor_norm > ci_pair_upper,
    significant_integrated = original_cor_int < ci_pair_integrated_lower | 
      original_cor_int > ci_pair_integrated_upper
  )
final_results
#write.csv(final_results, "output/NGA/boot_results_NGA.csv") #in the other script this ended up as v2


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
    labs(title = "Probability Density Function",
         subtitle = "Dashed line = observed correlation",
         x = "Correlation coefficient",
         y = "Density",
         caption = "NGA – N. cristatus – spring")
)

# export and save plot
#ggsave("figures/NGA/bootstrapPDF_NGA_Nc_spring.png", 
#       PDF_plot, 
#       width = 10, height = 6, dpi = 300, 
#       bg = "white")



## ------------------------------------------ ##
#     Plots -----
## ------------------------------------------ ##'

## ---------------- ##
#     Data -----
## ---------------- ##
PDO <- read.csv(file.path("output",
                          "NGA",
                          "PDO_NGA_integrated.csv"), 
                header = T) %>%
  select(-X) %>% 
  mutate(time = as.Date(time, 
                        format = "%Y-%m-%d"))

neocal <- read.csv(file.path("raw",
                             "NGA",
                             "NeocalanusCristBiomass_Spr.csv")) %>% #check version !!!!***
  mutate(Year = as.numeric(gsub("S", "", Year)),#remove S from year col
         taxa = "NeocalanusCrist",
         date = as.Date(paste0(Year, "-03-01"))) #give abundance daymonth to match w index

neocal$Yc_mean <- mean(neocal$LogMean)
neocal$Yc_sd <- sd(neocal$LogMean)
#zscore anomaly 
neocal$Anomaly_yr <- (neocal$LogMean - neocal$Yc_mean)/neocal$Yc_sd


cor_NGA <- read.csv(file.path("output",
                              "NGA",
                              "weightedCor_NGA.csv")) #check which version to use here

cor1 <- cor_NGA$cor_Norm[1]
cor2 <- cor_NGA$cor_Int[1]

pval1 <- final_results$pval_pair[1]
pval2 <- final_results$pval_pair_integrated[1]

(NormPlot <- ggplot() +
    # Driver time series (PDO normalized)
    geom_line(data = PDO, aes(x = time, y = pdoNorm, 
                              color = "Pacific Decadal Oscillation"), 
              linewidth = 1.2) +
    # Neocalanus cristatus time series (Biomass)
    geom_line(data = neocal, aes(x = date, y = Anomaly_yr, color = "Biomass"), 
              linewidth = 1.2) +
    # annotate correlation coefficient
    annotate("text", x = as.Date("2000-01-01"), y = 3,
             label = sprintf("r == %.2f", cor1),
             parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
    # annotate bootstrapped p-value
    annotate("text", x = as.Date("2000-01-01"), y = 2.4,
             label = if (pval1 < 0.001) {
               as.expression(bquote(italic(p) < 0.001))
             } else {
               as.expression(bquote(italic(p) == .(round(pval1, 3))))
             },
             parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
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
          plot.title = element_blank(),
          panel.grid = element_blank()) +
    # x-axis settings
    scale_x_date(breaks = seq(as.Date("1995-01-01"), as.Date("2022-01-01"), by = "5 years"),
                 date_labels = "%Y", expand = c(0, 0)) +
    coord_cartesian(xlim = as.Date(c("1995-01-01", "2022-01-01"))) +
    # manual color legend
    scale_color_manual(values = c("Pacific Decadal Oscillation" = "red", 
                                  "Biomass" = "blue")) + #confirm whether i have biomass or abundance
    labs(y = "PDO")
)

(Int1Plot <- ggplot() +
    # Driver time series (Integrated PDO)
    geom_line(data = PDO, aes(x = time, y = pdoInt, 
                              color = "Pacific Decadal Oscillation"), 
              linewidth = 1.2) +
    # Neocalanus cristatus time series (Biomass)
    geom_line(data = neocal, aes(x = date, y = Anomaly_yr, color = "Biomass"), 
              linewidth = 1.2) +
    # annotate correlation coefficient
    annotate("text", x = as.Date("2000-01-01"), y = 3,
             label = sprintf("r == %.2f", cor2),
             parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
    # annotate boostrapped p-value
    annotate("text", x = as.Date("2000-01-01"), y = 2.4,
             label = if (pval2 < 0.001) {
               as.expression(bquote(italic(p) < 0.001))
             } else {
               as.expression(bquote(italic(p) == .(round(pval2, 3))))
             },
             parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 19, color = "black", face = "bold"),
          axis.text = element_text(size = 15, color = "black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 15, face = "bold"),
          legend.key.size = unit(2, "lines"),
          plot.title = element_blank()) +
    # x-axis settings
    scale_x_date(breaks = seq(as.Date("1995-01-01"), as.Date("2022-01-01"), by = "5 years"),
                 date_labels = "%Y", expand = c(0, 0)) +
    coord_cartesian(xlim = as.Date(c("1995-01-01", "2022-01-01"))) +
    # manual color legend
    scale_color_manual(values = c("Pacific Decadal Oscillation" = "red", 
                                  "Biomass" = "blue")) +
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
  draw_label(label = bquote("NGA " * italic("Neocalanus cristatus") * " Spring"), 
             size = 18) +
  theme(plot.background = element_rect(fill = "transparent"))


# combine title and plots
(final_plot <- plot_grid(title, combined_plots, legend, ncol = 1, 
                         rel_heights = c(0.07, 0.8, 0.07))) 

# export and save plot
#ggsave("figures/NGA/boostrap_NGA_Nc_spring.png", 
#       final_plot, 
#       width = 8, height = 8, dpi = 300, 
#       bg = "white")
