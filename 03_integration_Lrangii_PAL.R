################################################################################
#############        LTER Pelagic Synthesis WG     #############################
#############   PAL - Double Integration Analysis  #############################
#############         Integration Calculations     #############################
## by: Alexandra Cabanelas 
## created FEB-2025, updated OCT-2025
################################################################################
## Palmer LTER
## Driver = Multivariate ENSO Index (MEI)
## Biology = Limacina rangii
## TAU = 1 year integration
## season = summer
## years = 1993 - 2023

# Script #3 : 03_integration_Lrangii_PAL

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
limacina <- read_csv(file.path("raw",
                               "PAL",
                               "PAL_Limacina.csv")) %>%
  rename(Anomaly_yr = Limacina) %>%
  mutate(taxa = "Limacina_rangii",
         date = as.Date(paste0(Year, "-03-01")), #give abundance daymonth to match w index
         # --- linear interpolation - these analyses dont like NAs -----
         #2021 & 2022 have NAs
         Anomaly_yr = approx(Year, Anomaly_yr, xout = Year)$y)
         
## ------------------------------------------ ##
#           2) Driver Data -----
## ------------------------------------------ ##
MEI <- read_csv(file.path("raw", 
                          "PAL", 
                          "MEI.csv")) %>% # contains up to FEB-2024
  mutate(
    #convert time to POSIXct
    time = as.POSIXct(as.Date(DATE, format = "%m/%d/%y")),
    mei = MEI,
    # --- Z-score MEI ----- 
    mei_z = (mei - mean(mei, na.rm = TRUE)) / sd(mei, na.rm = TRUE)
    #or mei_z = scale(mei)[, 1]
    ) %>%
  select(time, mei, mei_z)
  
mean(MEI$mei_z)  # should be close to 0
sd(MEI$mei_z)    # should be close to 1
hist(MEI$mei_z) 
plot(MEI$time, MEI$mei_z, type = "l")  

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
tau_bio <- 12        # biology memory in months
dt <- 1              # monthly time steps

## ------------------------------------------ ##
#  Apply integration (biological response) -----
## ------------------------------------------ ##

MEI <- MEI %>%
  mutate(
    mei_int = recursive_integration(mei_z, 
                                    tau = tau_bio, 
                                    dt = dt),
    # zscore after integration
    mei_int_z = (mei_int - mean(mei_int, na.rm = TRUE)) / sd(mei_int, na.rm = TRUE),  
    # or mei_int_z = scale(MEI$mei_int)[,1]  
    time = as.Date(time),
    tau_months = tau_bio)

## ------------------------------------------ ##
#   3) Test Correlations -----
## ------------------------------------------ ##

# --- interpolate (bio is yearly, mei is monthly) -----
limacina <- limacina %>%
  mutate(
    mei_orig = approx(MEI$time, MEI$mei_z, xout = date, rule = 1)$y,
    mei_int  = approx(MEI$time, MEI$mei_int_z, xout = date, rule = 1)$y
  )
#rule 2 = extend edge values beyond data range; rule 1 & 2 give same results here

# --- correlations -----
cor_raw <- cor(limacina$Anomaly_yr, limacina$mei_orig)
cor_int <- cor(limacina$Anomaly_yr, limacina$mei_int)

# --- get conventional p-values -----
pval_raw <- cor.test(limacina$Anomaly_yr, limacina$mei_orig)$p.value
pval_int <- cor.test(limacina$Anomaly_yr, limacina$mei_int)$p.value

## ------------------------------------------ ##
#   Store Results -----
## ------------------------------------------ ##
cor_PAL <- data.frame(site = "PAL",
                      taxa = "Limacina_rangii",
                      season = "summer",
                      tau_months = tau_bio,
                      cor_Norm = cor_raw, 
                      cor_Int = cor_int, 
                      pval_Norm = pval_raw, 
                      pval_Int = pval_int
)
cor_PAL

## ------------------------------------------ ##
#   Export CSVs -----
## ------------------------------------------ ##

# save L. rangii Z-score data 
#write.csv(limacina[, c("date", "Anomaly_yr", "taxa")], "output/PAL/Lrangii_anomalies.csv", row.names = FALSE)

# save MEI monthly z-scored
#write.csv(MEI, "output/PAL/MEI_int_Lrangii_PAL.csv", row.names = FALSE)

# save integration results 
#write.csv(cor_PAL, "output/PAL/cor_Lrangii_PAL.csv")

## ------------------------------------------ ##
#     Plots -----
## ------------------------------------------ ## 
(NormPlot <- ggplot() +
  # Driver time series (MEI normalized)
  geom_line(data = MEI, aes(x = time, y = mei_z, color = "MEI"), 
            linewidth = 1.2) +
  # Limacina time series (Abundance)
  geom_line(data = limacina, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
            linewidth = 1.2) +
  # annotate correlation coefficient
  annotate("text", x = as.Date("1993-01-01"), y = 3,
           label = sprintf("r == %.2f", cor_raw),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  # annotate p-value
  annotate("text", x = as.Date("1993-01-01"), y = 2.4,
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
  scale_x_date(breaks = seq(as.Date("1990-01-01"), as.Date("2025-01-01"), by = "10 years"),
               date_labels = "%Y", expand = c(0, 0)) +
  coord_cartesian(xlim = as.Date(c("1990-01-01", "2025-01-01"))) +
  # manual color legend
  scale_color_manual(values = c("MEI" = "red", "Abundance" = "blue")) +
  labs(y = "MEI")
)

(Int1Plot <- ggplot() +
  # Driver time series (Integrated MEI)
  geom_line(data = MEI, aes(x = time, y = mei_int_z, color = "MEI"), 
            linewidth = 1.2) +
  # Limacina time series (Abundance)
  geom_line(data = limacina, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
            linewidth = 1.2) +
  # annotate correlation coefficient
  annotate("text", x = as.Date("1993-01-01"), y = 3,
           label = sprintf("r == %.2f", cor_int),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  # annotate p-value
  annotate("text", x = as.Date("1993-01-01"), y = 2.4,
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
  scale_x_date(breaks = seq(as.Date("1990-01-01"), as.Date("2025-01-01"), by = "10 years"),
               date_labels = "%Y", expand = c(0, 0)) +
  coord_cartesian(xlim = as.Date(c("1990-01-01", "2025-01-01"))) +
  # manual color legend
  scale_color_manual(values = c("MEI" = "red", "Abundance" = "blue")) +
  # y-axis label
  labs(y = "Integrated MEI")
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
  draw_label(label = bquote("PAL " * italic("Limacina rangii") * " Summer"), 
             size = 18) +
  theme(plot.background = element_rect(fill = "transparent"))

# combine title and plots
(final_plot <- plot_grid(title, combined_plots, legend, ncol = 1, 
                        rel_heights = c(0.07, 0.8, 0.07))) 

# export and save plot
#ggsave("figures/PAL/conventional_PAL_Lr_summer.png", 
#       final_plot, 
#       width = 8, height = 8, dpi = 300, 
#       bg = "white")
##################################################################

library(cocor)
cocor_result <- cocor(~ Anomaly_yr + mei_orig | Anomaly_yr + mei_int, data = limacina_df)
print(cocor_result)


# correlation between the two predictors (critical for the test)
r_pred <- cor(limacina$mei_orig, limacina$mei_int,   use = "pairwise.complete.obs")
n      <- sum(complete.cases(limacina$Anomaly_yr, limacina$mei_orig, limacina$mei_int))

## --- Steiger/Williams test: are r_raw and r_int significantly different? ---
cc <- cocor.dep.groups.overlap(
  r.jk = cor_raw,         # cor(bio, driver_raw)
  r.jh = cor_int,         # cor(bio, driver_int)
  r.kh = r_pred,        # cor(driver_raw, driver_int)
  n = n,
  alternative = "two.sided"    # (reports multiple; steiger/williams are the standards)
)

## --- extract a simple summary ---
p_diff <- cc@steiger$`p.value`        # p-value for difference between the two r's
z_diff <- cc@steiger$`z`              # test statistic
r_diff <- r_int - r_raw

## --- quick R^2s ---
R2_raw <- r_raw^2
R2_int <- r_int^2
dR2    <- R2_int - R2_raw

## --- assemble tidy row (append to your cor_PAL) ---
cor_PAL_test <- tibble(
  site   = "PAL",
  taxa   = "Limacina_rangii",
  season = "summer",
  r_norm = r_raw,
  r_int  = r_int,
  R2_norm = R2_raw,
  R2_int  = R2_int,
  dR2     = dR2,
  r_pred  = r_pred,     # correlation between predictors
  n       = n,
  z_diff  = z_diff,
  p_diff  = p_diff
)

cor_PAL_test

cc <- cocor.dep.groups.overlap(
  r.jk = cor_raw,   # cor(bio, driver_raw)
  r.jh = cor_int,   # cor(bio, driver_int)
  r.kh = r_pred,    # cor(driver_raw, driver_int)
  n = n,
  alternative = "two.sided"
)
