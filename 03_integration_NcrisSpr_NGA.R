################################################################################
#############        LTER Pelagic Synthesis WG     #############################
#############   NGA - Double Integration Analysis  #############################
#############         Integration Calculations     #############################
## by: Alexandra Cabanelas 
## created FEB-2025, updated OCT-2025
################################################################################
## Northern Gulf of Alaska LTER
## Driver = Pacific Decadal Oscillation (PDO)
## Biology = N. cristatus
## TAU = 6 months integration time 
## season = spring 
## years = 1998 - 2022

# Script #3 : 03_integration_NcrisSpr_NGA

# script to integrate and test correlations

#######   ----    Double Integration Calculations:
# 1) Biological TS
#     a. Log(x + min/2) transformation
#     b. Spatial average
#     c. Z-score transformation for plotting on same axes
# 2) Driver TS 
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
neocal <- read_csv(file.path("raw",
                             "NGA",
                             "NeocalanusCristBiomass_Sprv2.csv")) %>%
  mutate(Year = as.numeric(gsub("S", "", Year)),#remove S from year col
         taxa = "NeocalanusCrist",
         date = as.Date(paste0(Year, "-03-01")), #give abundance daymonth to match w index
         # --- calculate anomalies -----
         Yc_mean = mean(LogMean, na.rm = TRUE),
         Yc_sd = sd(LogMean, na.rm = TRUE),
         Anomaly_yr = (LogMean - Yc_mean) / Yc_sd) 

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
tau_bio <- 6        # biology memory in months
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
neocal <- neocal %>%
  mutate(
    pdo_orig = approx(PDO$time, PDO$pdo_z, xout = date, rule = 1)$y,
    pdo_int  = approx(PDO$time, PDO$pdo_int_z, xout = date, rule = 1)$y
  )
#rule 2 = extend edge values beyond data range; rule 1 & 2 give same results here

# --- correlations -----
cor_raw <- cor(neocal$Anomaly_yr, neocal$pdo_orig)
cor_int <- cor(neocal$Anomaly_yr, neocal$pdo_int)

# --- get conventional p-values -----
pval_raw <- cor.test(neocal$Anomaly_yr, neocal$pdo_orig)$p.value
pval_int <- cor.test(neocal$Anomaly_yr, neocal$pdo_int)$p.value

## ------------------------------------------ ##
#   Store Results -----
## ------------------------------------------ ##

cor_NGA <- data.frame(site = "NGA",
                      taxa = "Neocalanus_cristatus",
                      season = "spring",
                      tau_months = tau_bio,
                      cor_Norm = cor_raw, 
                      cor_Int = cor_int, 
                      pval_Norm = pval_raw, 
                      pval_Int = pval_int
)
cor_NGA

## ------------------------------------------ ##
#   Export CSVs -----
## ------------------------------------------ ##

# save PDO monthly z-scored
#write.csv(PDO, "output/NGA/PDO_int_NcrisSpr_NGA.csv", row.names = FALSE)

# save integration results 
#write.csv(cor_NGA, "output/NGA/cor_NcrisSpr_NGA.csv")

## ------------------------------------------ ##
#     Plots -----
## ------------------------------------------ ## 
(NormPlot <- ggplot() +
  # Driver time series (PDO normalized)
  geom_line(data = PDO, aes(x = time, y = pdo_z, color = "PDO"), 
            linewidth = 1.2) +
  # Neocalanus cristatus time series (Abundance)
  geom_line(data = neocal, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
            linewidth = 1.2) +
  # annotate correlation coefficient
  annotate("text", x = as.Date("2000-01-01"), y = 3,
           label = sprintf("r == %.2f", cor_raw),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  # annotate p-value
  annotate("text", x = as.Date("2000-01-01"), y = 2.4,
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
  scale_x_date(breaks = seq(as.Date("1995-01-01"), as.Date("2022-01-01"), by = "5 years"),
               date_labels = "%Y", expand = c(0, 0)) +
  coord_cartesian(xlim = as.Date(c("1995-01-01", "2022-01-01"))) +
  # manual color legend
  scale_color_manual(values = c("PDO" = "red", "Abundance" = "blue")) +
  labs(y = "PDO")
)

(Int1Plot <- ggplot() +
  # Driver time series (Integrated PDO)
  geom_line(data = PDO, aes(x = time, y = pdo_int_z, color = "PDO"), 
            linewidth = 1.2) +
  # Neocalanus cristatus time series (Abundance)
  geom_line(data = neocal, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
            linewidth = 1.2) +
  # annotate correlation coefficient
  annotate("text", x = as.Date("2000-01-01"), y = 3,
           label = sprintf("r == %.2f", cor_int),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  # annotate p-value
  annotate("text", x = as.Date("2000-01-01"), y = 2.4,
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
  scale_x_date(breaks = seq(as.Date("1995-01-01"), as.Date("2022-01-01"), by = "5 years"),
               date_labels = "%Y", expand = c(0, 0)) +
  coord_cartesian(xlim = as.Date(c("1995-01-01", "2022-01-01"))) +
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
  draw_label(label = bquote("NGA " * italic("Neocalanus cristatus") * " Spring"), 
             size = 18) +
  theme(plot.background = element_rect(fill = "transparent"))

# combine title and plots
(final_plot <- plot_grid(title, combined_plots, legend, ncol = 1, 
                         rel_heights = c(0.07, 0.8, 0.07))) 

# export and save plot
#ggsave("figures/NGA/conventional_NGA_Nc_spring.png", 
#       final_plot, 
#       width = 8, height = 8, dpi = 300, 
#       bg = "white")
##################################################################