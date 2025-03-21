################################################################################
#############          Pelagic Synthesis           #############################
#############             FEB-2025                 #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas 
################################################################################
##### PAL LTER
##### Driver = Multivariate ENSO Index (MEI)
##### Biology = Limacina rangii
##### 1 year integration
##### season = winter
##### years = 1993 - 2023

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
#library(lubridate)
#library(zoo)

## ------------------------------------------ ##
#            Biology Data -----
## ------------------------------------------ ##
limacina <- read.csv(file.path("raw",
                               "PAL",
                               "PAL_Limacina.csv")) %>%
  mutate(date = as.Date(paste0(Year, "-03-01"))) #give abundance daymonth to match w index

#linear interpolation 
#these analyses dont like NAs - 2021 & 2022 have NAs
limacina$Limacina <- approx(limacina$Year, 
                            limacina$Limacina,
                            xout = limacina$Year)$y
#limacina$Limacina <- na.spline(limacina$Limacina)

## ------------------------------------------ ##
#            Driver Data -----
## ------------------------------------------ ##
MEI <- read.csv(file.path("raw", 
                          "PAL", 
                          "MEI.csv")) %>%
  mutate(time = as.POSIXct(as.Date(DATE, format = "%m/%d/%y")),
         mei = MEI) %>%
  select(time, mei) %>%
  as.data.frame()

#write.csv(MEI, "output/PAL/MEI_PAL_integrated.csv")
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

calculateIntegrations = function(data, tau = 365*1, 
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
    
    par(mfrow = c(2,1))
    plot(data[,1], data[[n]], type = 'l', ylab = n, xlab = '') #original signal
    plot(data[,1], data[[paste0(n,'Int')]], type = 'l', ylab = paste0(n, 'Int'), xlab = '') #integration1
    #plot(data[,1], data[[paste0(n, 'DInt')]], type = 'l', ylab = paste0(n, ' DInt')) #integration2
  }
  data
}

meiInt <- calculateIntegrations(MEI) #if your R window is too small may give err
meiInt$time <- as.Date(meiInt$time)

## ------------------------------------------ ##
#     Correlations -----
## ------------------------------------------ ##
# 3) Test correlations between bio TS + each driver TS 

driver1 <- meiInt[, c(1,3)]
int1 <- meiInt[, c(1,4)]
#int2 <- meiInt[, c(1,5)]

# interpolate driver values - matching dates/times 
interpolated_driver <- approx(driver1$time, driver1[[2]], xout = limacina$date)$y
interpolated_int1 <- approx(int1$time, int1[[2]], xout = limacina$date)$y
#interpolated_int2 <- approx(int2$time, int2[[2]], xout = limacina$date)$y

# calculate correlations
cor1 <- cor(limacina$Limacina, interpolated_driver)
cor2 <- cor(limacina$Limacina, interpolated_int1)
#cor3 <- cor(limacina$Limacina, interpolated_int2)

# get conventional pvalues 
pval1 <- cor.test(limacina$Limacina, interpolated_driver)$p.value
pval2 <- cor.test(limacina$Limacina, interpolated_int1)$p.value
#pval3 <- cor.test(limacina$Limacina, interpolated_int2)$p.value

wcor_PAL <- data.frame(site = "PAL",
                       taxa = "limacina",
                       season = "winter",
                       cor_Norm = cor1, 
                       cor_Int = cor2, 
                       #cor_DouInt = cor_DouInt, 
                       pval_Norm = pval1, 
                       pval_Int = pval2
                       #pval_DouInt = pval_DouInt
)
wcor_PAL
#write.csv(wcor_PAL, "output/PAL/weightedCor_PAL.csv")
## ------------------------------------------ ##
#     Plots -----
## ------------------------------------------ ## 
(NormPlot <- ggplot() +
  # Driver time series (MEI normalized)
  geom_line(data = driver1, aes(x = time, y = meiNorm, color = "MEI"), 
            linewidth = 1.2) +
  # Limacina time series (Abundance)
  geom_line(data = limacina, aes(x = date, y = Limacina, color = "Abundance"), 
            linewidth = 1.2) +
  # annotate correlation coefficient
  annotate("text", x = as.Date("1990-01-01"), y = 3,
           label = sprintf("rho == %.4f", cor1),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  # annotate p-value
  annotate("text", x = as.Date("1990-01-01"), y = 2.4,
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
  scale_x_date(breaks = seq(as.Date("1990-01-01"), as.Date("2025-01-01"), by = "10 years"),
               date_labels = "%Y", expand = c(0, 0)) +
  coord_cartesian(xlim = as.Date(c("1990-01-01", "2025-01-01"))) +
  # manual color legend
  scale_color_manual(values = c("MEI" = "red", "Abundance" = "blue")) +
  labs(y = "MEI")
)

(Int1Plot <- ggplot() +
  # Driver time series (Integrated MEI)
  geom_line(data = int1, aes(x = time, y = meiInt, color = "MEI"), 
            linewidth = 1.2) +
  # Limacina time series (Abundance)
  geom_line(data = limacina, aes(x = date, y = Limacina, color = "Abundance"), 
            linewidth = 1.2) +
  # annotate correlation coefficient
  annotate("text", x = as.Date("1990-01-01"), y = 3,
           label = sprintf("rho == %.4f", cor2),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  # annotate p-value
  annotate("text", x = as.Date("1990-01-01"), y = 2.4,
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
ggsave("figures/PAL/conventional_PAL_Lr_summer.png", 
       final_plot, 
       width = 8, height = 8, dpi = 300, 
       bg = "white")
##################################################################


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
  scale_x_date(breaks = seq(as.Date("1990-01-01"), as.Date("2025-01-01"), by = "10 years"),
               date_labels = "%Y", expand = c(0, 0)) +
  coord_cartesian(xlim = as.Date(c("1990-01-01", "2025-01-01"))) +
  
  # Manual color legend
  scale_color_manual(values = c("MEI" = "red", "Abundance" = "blue")) +
  
  # Y-axis label
  labs(y = "Double Integrated MEI")

Int2Plot <- Int2Plot + theme(legend.position = "none")
