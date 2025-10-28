################################################################################
#############        LTER Pelagic Synthesis WG     #############################
#############   CCE - Double Integration Analysis  #############################
#############         Integration Calculations     #############################
## by: Alexandra Cabanelas 
## created FEB-2025, updated OCT-2025
################################################################################
## California Current Ecosystem LTER
## Driver = Pacific Decadal Oscillation (PDO)
## Biology = Nyctiphanes simplex (abundance m2)
## TAU = 2-year integration time 
## season = spring (Feb, Mar, Apr, May)
## years = 1951 - 2021
## night tows only

# Script #3 : 03_integration_Nsimplex_CCE

# script to integrate and test correlations

#######   ----   Integration Calculations:
# 1) Biological TS (annual)
#     a. Log(x + min/2) transformation == for CCE, Log(Abundance per m2 + 1)
#     b. Spatial average
#     c. Z-score transformation for plotting on same axes
# 2) Driver TS (monthly)
#     a. Raw TS (at higher resolution than bio TS)
#     b. Temporally average (backwards in time) over characteristic 
#        time span of organism - 1st integration
#     c. Temporally average again - 2nd integration 
# INTEGRATING ONCE IN THIS ANALYSIS*
# 3) Test correlations between bio TS + each driver TS 
################################################################################

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(tidyverse) #v2.0.0
library(cowplot) #v1.1.3

## ------------------------------------------ ##
#           1) Biology Data -----
## ------------------------------------------ ##
euphs <- read_csv(file.path("raw",
                            "CCE",
                            "nsimplex_CCE.csv")) %>%
  mutate(taxa = "Nsimplex") %>% 
  # --- add the missing years -----
  #no sampling = 1967, 1968, 1971, 1973 and 2020
  #real 0s = 1972, 1976, 2010-2012
  complete(Year = seq(min(Year), max(Year), by = 1)) %>%
  arrange(Year) %>%
  mutate(
    # --- linear interpolation - these analyses dont like NAs -----
    Abundance = approx(Year, Abundance, xout = Year)$y,  
  
    # --- calculate anomalies -----
    # Abundance = Log10(Abundance per m2 + 1)
    Yc_mean = mean(Abundance, na.rm = TRUE),
    # the matlab code uses sample (denominator = n-1) sd
    Yc_sd = sd(Abundance, na.rm = TRUE), 
    #^same as: sqrt(sum((Abundance - Yc_mean)^2, na.rm = TRUE) / (sum(!is.na(Abundance)) - 1))
    #for pop sd(denominator=n) = sqrt(sum((Abundance - Yc_mean)^2, na.rm = TRUE) / sum(!is.na(Abundance)))
  
    # --- Z-score ----- 
    Anomaly_yr = (Abundance - Yc_mean) / Yc_sd, 
    # or Anomaly_yr = scale(Abundance)[, 1] #sample SD
    
    # --- give abundance data day & month to match index (March 1st to each yr)
    date = as.Date(paste0(Year, "-03-01"))                
  ) %>%
  # --- fill region and taxa metadata so they arent NA -----
  fill(Region, taxa, .direction = "downup")              

## ------------------------------------------ ##
#           2) Driver Data -----
## ------------------------------------------ ##
PDO <- read_csv(file.path("raw",
                          "CCE",
                          "PDO.csv")) %>% # contains up to DEC-2024
  # long format
  pivot_longer(cols = Jan:Dec,
               names_to = "month",
               values_to = "pdo") %>%
  #filtering 10 yrs before bio data
  filter(Year >= 1941, Year <= 2022, pdo < 99) %>%
  mutate(
    #turn month into integer
    month = match(month, month.abb),
    #convert time to POSIXct
    time = as.POSIXct(paste(Year, month, "01"), 
                      format = "%Y %m %d", 
                      tz = "UTC"),
    # --- Z-score PDO ----- 
    pdo_z = (pdo - mean(pdo, na.rm = TRUE)) / sd(pdo, na.rm = TRUE)
    #or pdo_z = scale(pdo)[, 1]
  ) %>%
  #as.data.frame() %>%
  select(time, pdo, pdo_z)

mean(PDO$pdo_z)  # should be close to 0
sd(PDO$pdo_z)    # should be close to 1
hist(PDO$pdo_z) 
plot(PDO$time, PDO$pdo_z, type = "l")

## ------------------------------------------ ##
#  Function for integrations (matching MATLAB/Manu's approach) -----
## ------------------------------------------ ##
# AR(1) style autoregressive smoothing
recursive_integration <- function(x, tau, dt = 1) {
  # calculate autoregressive weight (how much past value is retained)
  # tau = memory in months, dt = time step 1 month
  alpha <- 1 - dt / tau # autoregressive weight; alpha bio, Euler
  y <- numeric(length(x)) #bio response
  y[1] <- x[1] # initial condition; (MATLAB: o(1) = forc.sig(1))
  for (i in 1:(length(x) - 1)) {
    y[i + 1] <- alpha * y[i] + x[i] * dt #Euler forward recursion loop
  }
  return(y)
}

# set parameters
tau_bio <- 24        # biology memory in months
dt <- 1              # monthly time steps

## ------------------------------------------ ##
#  Apply integration (biological response) -----
## ------------------------------------------ ##

PDO <- PDO %>%
  mutate(
    pdo_int = recursive_integration(pdo_z, 
                                    tau = tau_bio, 
                                    dt = dt),
    # zscore after integration
    pdo_int_z = (pdo_int - mean(pdo_int, na.rm = TRUE)) / sd(pdo_int, na.rm = TRUE),  
    # or pdo_int_z = scale(PDO$pdo_int)[,1]  
    time = as.Date(time),
    tau_months = tau_bio)

## ------------------------------------------ ##
#   3) Test Correlations -----
## ------------------------------------------ ##

# --- interpolate (bio is yearly, pdo is monthly) -----
euphs <- euphs %>%
  mutate(
    pdo_orig = approx(PDO$time, PDO$pdo_z, xout = date, rule = 1)$y,
    pdo_int  = approx(PDO$time, PDO$pdo_int_z, xout = date, rule = 1)$y
  )
#rule 2 = extend edge values beyond data range; rule 1 & 2 give same results here

# --- correlations -----
cor_raw <- cor(euphs$Anomaly_yr, euphs$pdo_orig)
cor_int <- cor(euphs$Anomaly_yr, euphs$pdo_int)

# --- get conventional p-values -----
pval_raw <- cor.test(euphs$Anomaly_yr, euphs$pdo_orig)$p.value
pval_int <- cor.test(euphs$Anomaly_yr, euphs$pdo_int)$p.value

## ------------------------------------------ ##
#   Store Results -----
## ------------------------------------------ ##
cor_CCE <- data.frame(site = "CCE",
                      taxa = "Nsimplex",
                      season = "spring",
                      tau_months = tau_bio,
                      cor_Norm = cor_raw, 
                      cor_Int = cor_int, 
                      pval_Norm = pval_raw, 
                      pval_Int = pval_int
)
cor_CCE

## ------------------------------------------ ##
#   Export CSVs -----
## ------------------------------------------ ##

# save N. simplex Z-score data 
#write.csv(euphs[, c("date", "Anomaly_yr", "taxa")], "output/CCE/Nsimplex_anomalies.csv", row.names = FALSE)

# save PDO monthly z-scored
#write.csv(PDO, "output/CCE/PDO_int_Nsimplex_CCE.csv", row.names = FALSE)

# save integration results 
#write.csv(cor_CCE, "output/CCE/cor_Nsimplex_CCE.csv")

## ------------------------------------------ ##
#     Plots -----
## ------------------------------------------ ## 
(NormPlot <- ggplot() +
  # Driver time series (PDO normalized)
  geom_line(data = PDO, aes(x = time, y = pdo_z, color = "PDO"), 
            linewidth = 1.2) +
  # Nsimplex time series (Abundance)
  geom_line(data = euphs, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
            linewidth = 1.2) +
  # annotate correlation coefficient
  annotate("text", x = as.Date("1960-01-01"), y = 3,
           label = sprintf("r == %.2f", cor_raw),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  # annotate p-value
  annotate("text", x = as.Date("1960-01-01"), y = 2.4,
           label = if (pval_raw < 0.001) {
             as.expression(bquote(italic(p) < 0.001))
           } else {
             as.expression(bquote(italic(p) == .(round(pval_raw, 3))))
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
  # annotate p-value
  annotate("text", x = as.Date("1960-01-01"), y = 2.4,
           label = if (pval_int < 0.001) {
             as.expression(bquote(italic(p) < 0.001))
           } else {
             as.expression(bquote(italic(p) == .(round(pval_int, 3))))
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
#ggsave("figures/CCE/conventional_CCE_Ns_spring.png", 
#       final_plot, 
#       width = 8, height = 8, dpi = 300, 
#       bg = "white")

plot(euphs$date, euphs$Anomaly_yr, type = "l", col = "black", lwd = 2)
lines(euphs$date, euphs$pdo_orig, col = "blue", lwd = 2)
lines(euphs$date, euphs$pdo_int, col = "red", lwd = 2)
legend("topright", legend = c("Biological anomaly", "PDO (raw)", "PDO (integrated)"),
       col = c("black", "blue", "red"), lwd = 2)
