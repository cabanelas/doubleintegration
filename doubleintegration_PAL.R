################################################################################
#############          Pelagic Synthesis           #############################
#############             FEB-2025                 #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas 
################################################################################
##### PAL LTER
##### Driver = Multivariate ENSO Index (MEI)
##### Biology = Limacina
## 1 year integration

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
# 3) Test correlations between bio TS + each driver TS 
################################################################################

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##

library(tidyverse) #v2.0.0
library(gridExtra) #v2.3 
library(cowplot) #v1.1.3
library(here) #v1.0.1; easily build path to files  
library(listviewer) #v4.0.0; for looking at lists interactively 
library(lubridate)
library(zoo)

## ------------------------------------------ ##
#            Biology Data -----
## ------------------------------------------ ##

limacina <- read.csv(file.path("raw",
                               "PAL",
                               "PAL_Limacina.csv"))

## ------------------------------------------ ##
#            Driver Data -----
## ------------------------------------------ ##
MEI <- read.csv(file.path("raw",
                          "PAL",
                          "MEI.csv"))

MEIa <- MEI %>%
  mutate(time = as.Date(DATE, format = "%m/%d/%y")) %>%
  select(time, MEI) %>%
  rename(mei = MEI)

# convert time column to POSIXct format
MEIa <- MEIa %>% 
  mutate(time = as.POSIXct(time, 
                           format = "%Y-%m-%d"))

MEIa <- as.data.frame(MEIa)
## ------------------------------------------ ##
#     Function for integrations -----
## ------------------------------------------ ##
# 2) Driver TS 
#     a. Raw TS (at higher resolution than bio TS)
#     b. Temporally average (backwards in time) over characteristic 
#        time span of organism - 1st integration
#     c. Temporally average again - 2nd integration 

calculateIntegrations = function(data, tau = 365*1, f = function(x) {mean(x, na.rm = T)}) {
  for (n in names(data)[-1]) {
    data[[paste0(n,'Norm')]] = (data[[n]] - mean(data[[n]], na.rm = T)) / sd(data[[n]]) #norm
    
    data[[paste0(n,'Int')]] = NA
    data[[paste0(n,'DInt')]] = NA
    
    for (i in 2:nrow(data)) {
      k = data[,1] <= data[i,1] & data[,1] > data[i,1] - tau * 86400 #secs
      data[[paste0(n,'Int')]][i] = f(data[[paste0(n,'Norm')]][k])
      data[[paste0(n,'DInt')]][i] = f(data[[paste0(n,'Int')]][k])
    }
    
    par(mfrow = c(3,1))
    plot(data[,1], data[[n]], type = 'l', ylab = n, xlab = '') #original signal
    plot(data[,1], data[[paste0(n,'Int')]], type = 'l', ylab = paste0(n, 'Int'), xlab = '') #integration1
    plot(data[,1], data[[paste0(n,'DInt')]], type = 'l', ylab = paste0(n, 'DInt')) #integration2
  }
  data
}

meiInt <- calculateIntegrations(MEIa)
meiInt$time <- as.Date(meiInt$time)

## ------------------------------------------ ##
#     Correlations -----
## ------------------------------------------ ##
# 3) Test correlations between bio TS + each driver TS 

driver1 <- meiInt[, c(1,3)]
int1 <- meiInt[, c(1,4)]
int2 <- meiInt[, c(1,5)]

# need to give abundance data random day & month just to match things
# DOUBLE CHECK THIS

limacina$date <- as.Date(paste0(limacina$Year, "-03-01"))


# Interpolate driver values - need to check it makes sense to use abundance data date 
interpolated_driver <- approx(driver1$time, driver1[[2]], xout = limacina$date)$y
interpolated_int1 <- approx(int1$time, int1[[2]], xout = limacina$date)$y
interpolated_int2 <- approx(int2$time, int2[[2]], xout = limacina$date)$y

limacina$Limacina <- na.spline(limacina$Limacina)

cor1 <- cor(limacina$Limacina, interpolated_driver)
cor2 <- cor(limacina$Limacina, interpolated_int1)
cor3 <- cor(limacina$Limacina, interpolated_int2)

pval1 <- cor.test(limacina$Limacina, interpolated_driver)$p.value
pval2 <- cor.test(limacina$Limacina, interpolated_int1)$p.value
pval3 <- cor.test(limacina$Limacina, interpolated_int2)$p.value



## ------------------------------------------ ##
#     Plots -----
## ------------------------------------------ ## 

NormPlot <- ggplot() +
  # Driver time series (MEI normalized)
  geom_line(data = driver1, aes(x = time, y = meiNorm, color = "MEI"), 
            linewidth = 1.2) +
  
  # LOESS smoothing for driver data
  #geom_smooth(data = driver1, aes(x = time, y = meiNorm), 
  #            method = "loess", span = 0.1, 
  #            se = FALSE, linewidth = 1, color = "black") +
  
  # Limacina time series (Abundance)
  geom_line(data = limacina, aes(x = date, y = Limacina, color = "Abundance"), 
            linewidth = 1.2) +
  
  # Annotate correlation coefficient
  annotate("text", x = as.Date("1990-01-01"), y = 3,
           label = sprintf("rho == %.4f", cor1),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  
  # Annotate p-value
  annotate("text", x = as.Date("1990-01-01"), y = 2.4,
           label = if (pval1 < 0.001) {
             as.expression(bquote(italic(P) < 0.001))
           } else {
             as.expression(bquote(italic(P) == .(round(pval1, 3))))
           },
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  
  # Theme settings
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
  
  # X-axis settings
  scale_x_date(breaks = seq(as.Date("1978-01-01"), as.Date("2025-01-01"), by = "10 years"),
               date_labels = "%Y", expand = c(0, 0)) +
  coord_cartesian(xlim = as.Date(c("1978-01-01", "2025-01-01"))) +
  
  # Manual color legend
  scale_color_manual(values = c("MEI" = "red", "Abundance" = "blue")) +
  labs(y = "MEI")


Int1Plot <- ggplot() +
  # Driver time series (Integrated MEI)
  geom_line(data = int1, aes(x = time, y = meiInt, color = "MEI"), 
            linewidth = 1.2) +
  
  # Limacina time series (Abundance)
  geom_line(data = limacina, aes(x = date, y = Limacina, color = "Abundance"), 
            linewidth = 1.2) +
  
  # Annotate correlation coefficient
  annotate("text", x = as.Date("1990-01-01"), y = 3,
           label = sprintf("rho == %.4f", cor2),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  
  # Annotate p-value
  annotate("text", x = as.Date("1990-01-01"), y = 2.4,
           label = if (pval2 < 0.001) {
             as.expression(bquote(italic(P) < 0.001))
           } else {
             as.expression(bquote(italic(P) == .(round(pval2, 3))))
           },
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  
  # Theme settings
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 19, color = "black", face = "bold"),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, face = "bold"),
        legend.key.size = unit(2, "lines"),
        plot.title = element_blank()) +
  
  # X-axis settings
  scale_x_date(breaks = seq(as.Date("1978-01-01"), as.Date("2025-01-01"), by = "10 years"),
               date_labels = "%Y", expand = c(0, 0)) +
  coord_cartesian(xlim = as.Date(c("1978-01-01", "2025-01-01"))) +
  
  # Manual color legend
  scale_color_manual(values = c("MEI" = "red", "Abundance" = "blue")) +
  
  # Y-axis label
  labs(y = "Integrated MEI")





Int2Plot <- ggplot() +
  # Driver time series (Integrated MEI)
  geom_line(data = int2, aes(x = time, y = meiDInt, color = "MEI"), 
            linewidth = 1.2) +
  
  # Limacina time series (Abundance)
  geom_line(data = limacina, aes(x = date, y = Limacina, color = "Abundance"), 
            linewidth = 1.2) +
  
  # Annotate correlation coefficient
  annotate("text", x = as.Date("1990-01-01"), y = 3,
           label = sprintf("rho == %.4f", cor3),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  
  # Annotate p-value
  annotate("text", x = as.Date("1990-01-01"), y = 2.4,
           label = if (pval3 < 0.001) {
             as.expression(bquote(italic(P) < 0.001))
           } else {
             as.expression(bquote(italic(P) == .(round(pval3, 3))))
           },
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  
  # Theme settings
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 19, color = "black", face = "bold"),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, face = "bold"),
        legend.key.size = unit(2, "lines"),
        plot.title = element_blank()) +
  
  # X-axis settings
  scale_x_date(breaks = seq(as.Date("1978-01-01"), as.Date("2025-01-01"), by = "10 years"),
               date_labels = "%Y", expand = c(0, 0)) +
  coord_cartesian(xlim = as.Date(c("1978-01-01", "2025-01-01"))) +
  
  # Manual color legend
  scale_color_manual(values = c("MEI" = "red", "Abundance" = "blue")) +
  
  # Y-axis label
  labs(y = "Double Integrated MEI")




legend <- cowplot::get_plot_component(NormPlot, 'guide-box-bottom', 
                                      return_all = TRUE)
cowplot::ggdraw(legend)

# Remove legend
NormPlot <- NormPlot + theme(legend.position = "none")
Int1Plot <- Int1Plot + theme(legend.position = "none")
Int2Plot <- Int2Plot + theme(legend.position = "none")

# Combine plots into one grid
combined_plots <- plot_grid(NormPlot, Int1Plot, Int2Plot,
                            ncol = 1, align = "v")

# Add plot title 
title <- ggdraw() +
  draw_label(label = bquote("PAL " * italic("Limacina rangii") * " Winter"), 
             size = 18) +
  theme(plot.background = element_rect(fill = "transparent"))


# Combine title and plots
final_plot <- plot_grid(title, combined_plots, legend, ncol = 1, 
                        rel_heights = c(0.07, 0.8, 0.07)) 
