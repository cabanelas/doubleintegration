################################################################################
#############          Pelagic Synthesis           #############################
#############             FEB-2025                 #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas 
################################################################################
##### NGA LTER
##### Driver = Pacific Decadal Oscillation (PDO)
##### Biology = Neocalanus cristatus
##### 6 months? integration time 
##### season = spring 
##### years = 1998 - 2022

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
#library(gridExtra) #v2.3 
library(cowplot) #v1.1.3
#library(here) #v1.0.1; easily build path to files  
#library(listviewer) #v4.0.0; for looking at lists interactively 
#library(lubridate)
#library(zoo)

## ------------------------------------------ ##
#            Biology Data -----
## ------------------------------------------ ##
neocal <- read.csv(file.path("raw",
                            "NGA",
                            "NeocalanusCristBiomass_Spr.csv")) %>%
  mutate(Year = as.numeric(gsub("S", "", Year)),#remove S from year col
         #Year = as.numeric(gsub("[A-Za-z]", "", Year))
         taxa = "NeocalanusCrist",
         date = as.Date(paste0(Year, "-03-01"))) #give abundance daymonth to match w index

## ------------------------------------------ ##
#            Tidy -----
## ------------------------------------------ ##
neocal$Yc_mean <- mean(neocal$LogMean)
neocal$Yc_sd <- sd(neocal$LogMean)
#zscore anomaly 
neocal$Anomaly_yr <- (neocal$LogMean - neocal$Yc_mean)/neocal$Yc_sd

## ------------------------------------------ ##
#            Driver Data -----
## ------------------------------------------ ##
PDO <- read.csv(file.path("raw",
                          "CCE",
                          "PDO.csv")) %>%
  pivot_longer(cols = Jan:Dec,
               names_to = "month",
               values_to = "pdo") %>%
  filter(pdo < 99 & Year > 1940) %>%
  mutate(
    month = match(month, month.abb),
    time = as.POSIXct(paste(Year, month, "01"), format = "%Y %m %d", tz = "UTC")  # Convert to POSIXct
  ) %>%
  as.data.frame() %>%
  select(time, pdo)

## ------------------------------------------ ##
#     Function for weighted integrations -----
## ------------------------------------------ ##
# 2) Driver TS 
#     a. Raw TS (at higher resolution than bio TS)
#     b. Temporally average (backwards in time) over characteristic 
#        time span of organism - 1st integration
#     c. Temporally average again - 2nd integration 

# weighted average
# exponential weighing - weights decrease for older data points
# w = e^(t-t_i)/lamba; where lambda controls the decay

calculateIntegrations = function(data, tau = 180, 
                                 f = function(x, w) {sum(x * w, na.rm = T) / sum(w, na.rm = T)}) {
  for (n in names(data)[-1]) {
    data[[paste0(n,'Norm')]] = (data[[n]] - mean(data[[n]], na.rm = T)) / sd(data[[n]]) # norm
    
    data[[paste0(n,'Int')]] = NA
    data[[paste0(n,'DInt')]] = NA
    
    for (i in 2:nrow(data)) {
      #first integration
      k1 = data[,1] <= data[i,1] & data[,1] > data[i,1] - tau * 86400 # secs
      tdiff1 = as.numeric(data[i,1] - data[k1,1]) # convert difftime to numeric
      w1 = exp(-tdiff1 / (tau * 86400 / log(2))) # weights
      data[[paste0(n,'Int')]][i] = f(data[[paste0(n,'Norm')]][k1], w1)
      
      # second integration
      k2 = data[,1] <= data[i,1] & data[,1] > data[i,1] - tau * 86400
      tdiff2 = as.numeric(data[i,1] - data[k2,1])
      w2 = exp(-tdiff2 / (tau * 86400 / log(2)))
      data[[paste0(n, 'DInt')]][i] = f(data[[paste0(n, 'Int')]][k2], w2)
    }
    
    par(mfrow = c(3,1))
    plot(data[,1], data[[n]], type = 'l', ylab = n, xlab = '') #original signal
    plot(data[,1], data[[paste0(n,'Int')]], type = 'l', ylab = paste0(n, 'Int'), xlab = '') #integration1
    plot(data[,1], data[[paste0(n, 'DInt')]], type = 'l', ylab = paste0(n, ' DInt')) #integration2
  }
  data
}

pdoInt <- calculateIntegrations(PDO) #if your R window is too small may give err
pdoInt$time <- as.Date(pdoInt$time)
#write.csv(pdoInt, "output/NGA/PDO_NGA_integrated.csv")
## ------------------------------------------ ##
#     Correlations -----
## ------------------------------------------ ##
# 3) Test correlations between bio TS + each driver TS 

driver1 <- pdoInt[, c(1,3)]
int1 <- pdoInt[, c(1,4)]
#int2 <- pdoInt[, c(1,5)]

# interpolate driver values - matching dates/times 
interpolated_driver <- approx(driver1$time, driver1[[2]], xout = neocal$date)$y
interpolated_int1 <- approx(int1$time, int1[[2]], xout = neocal$date)$y
#interpolated_int2 <- approx(int2$time, int2[[2]], xout = neocal$date)$y

# calculate correlations
cor1 <- cor(neocal$Anomaly_yr, interpolated_driver)
cor2 <- cor(neocal$Anomaly_yr, interpolated_int1)
#cor3 <- cor(neocal$Anomaly_yr, interpolated_int2)

# get conventional pvalues
pval1 <- cor.test(neocal$Anomaly_yr, interpolated_driver)$p.value
pval2 <- cor.test(neocal$Anomaly_yr, interpolated_int1)$p.value
#pval3 <- cor.test(neocal$Anomaly_yr, interpolated_int2)$p.value

wcor_NGA <- data.frame(site = "NGA",
                       taxa = "Neocalanus_cristatus",
                       season = "spring",
                       cor_Norm = cor1, 
                       cor_Int = cor2, 
                       #cor_DouInt = cor_DouInt, 
                       pval_Norm = pval1, 
                       pval_Int = pval2
                       #pval_DouInt = pval_DouInt
)
wcor_NGA
#write.csv(wcor_NGA, "output/NGA/weightedCor_NGA.csv")

## ------------------------------------------ ##
#     Plots -----
## ------------------------------------------ ## 
(NormPlot <- ggplot() +
  # Driver time series (PDO normalized)
  geom_line(data = driver1, aes(x = time, y = pdoNorm, color = "PDO"), 
            linewidth = 1.2) +
  # Neocalanus cristatus time series (Abundance)
  geom_line(data = neocal, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
            linewidth = 1.2) +
  # annotate correlation coefficient
  annotate("text", x = as.Date("2000-01-01"), y = 3,
           label = sprintf("r == %.2f", cor1),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  # annotate p-value
  annotate("text", x = as.Date("2000-01-01"), y = 2.4,
           label = if (pval1 < 0.001) {
             as.expression(bquote(italic(P) < 0.001))
           } else {
             as.expression(bquote(italic(P) == .(round(pval1, 3))))
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
  geom_line(data = int1, aes(x = time, y = pdoInt, color = "PDO"), 
            linewidth = 1.2) +
  # Neocalanus cristatus time series (Abundance)
  geom_line(data = neocal, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
            linewidth = 1.2) +
  # annotate correlation coefficient
  annotate("text", x = as.Date("2000-01-01"), y = 3,
           label = sprintf("r == %.2f", cor2),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  # annotate p-value
  annotate("text", x = as.Date("2000-01-01"), y = 2.4,
           label = if (pval2 < 0.001) {
             as.expression(bquote(italic(P) < 0.001))
           } else {
             as.expression(bquote(italic(P) == .(round(pval2, 3))))
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
ggsave("figures/NGA/conventional_NGA_Nc_spring.png", 
       final_plot, 
       width = 8, height = 8, dpi = 300, 
       bg = "white")
##################################################################


Int2Plot <- ggplot() +
  # Driver time series (Integrated PDO)
  geom_line(data = int2, aes(x = time, y = pdoDInt, color = "MEI"), 
            linewidth = 1.2) +
  # Neocalanus cristatus time series (Abundance)
  geom_line(data = neocal, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
            linewidth = 1.2) +
  # annotate correlation coefficient
  annotate("text", x = as.Date("2000-01-01"), y = 3,
           label = sprintf("r == %.2f", cor3),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  # annotate p-value
  annotate("text", x = as.Date("2000-01-01"), y = 2.4,
           label = if (pval3 < 0.001) {
             as.expression(bquote(italic(P) < 0.001))
           } else {
             as.expression(bquote(italic(P) == .(round(pval3, 3))))
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
  scale_x_date(breaks = seq(as.Date("1995-01-01"), as.Date("2022-01-01"), by = "10 years"),
               date_labels = "%Y", expand = c(0, 0)) +
  coord_cartesian(xlim = as.Date(c("1995-01-01", "2022-01-01"))) +
  # manual color legend
  scale_color_manual(values = c("PDO" = "red", "Abundance" = "blue")) +
  # y-axis label
  labs(y = "Double Integrated PDO")

Int2Plot <- Int2Plot + theme(legend.position = "none")