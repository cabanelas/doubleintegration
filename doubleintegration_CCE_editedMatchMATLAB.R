################################################################################
#############          Pelagic Synthesis           #############################
#############             FEB-2025                 #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas 
################################################################################
##### CCE LTER
##### Driver = Pacific Decadal Oscillation (PDO)
##### Biology = Nyctiphanes simplex
##### 2 year integration time 
##### season = spring 
##### years = 1951 - 2021

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
library(cowplot) #v1.1.3

## ------------------------------------------ ##
#            Biology Data -----
## ------------------------------------------ ##
euphs <- read.csv(file.path("raw",
                            "CCE",
                            "nsimplex_CCE.csv")) %>%
  mutate(taxa = "Nsimplex") 

## ------------------------------------------ ##
#            Interpolation -----
## ------------------------------------------ ##
#need to add the missing years 
#no sampling 1967, 1968, 1971, 1973 and 2020
#real 0s = 1972, 1976, 2010-2012
full_years <- data.frame(Year = seq(min(euphs$Year), max(euphs$Year), by = 1))

euphs <- full_years %>%
  left_join(euphs, by = "Year")

rm(full_years)

#linear interpolation 
#these analyses dont like NAs
euphs$Abundance <- approx(euphs$Year, 
                          euphs$Abundance, 
                          xout = euphs$Year)$y

## ------------------------------------------ ##
#            Zscore -----
## ------------------------------------------ ##
# the matlab code does use sample sd
meanNaN <- function(x) mean(x, na.rm = TRUE)
sdNaN <- function(x) sd(x, na.rm = TRUE)  # sample SD

#euphs$Yc_mean <- mean(euphs$Abundance)
#euphs$Yc_sd <- sqrt(sum((euphs$Abundance - mean(euphs$Abundance, na.rm = TRUE))^2, na.rm = TRUE) / 
#       sum(!is.na(euphs$Abundance))) #pop mean; not sample..


#another option is to use sample SD by doing
# euphs$Anomaly_yr <- scale(euphs$Abundance)[,1]  # uses sample SD (n-1), not population

#zscored
#euphs$Anomaly_yr <- (euphs$Abundance - euphs$Yc_mean)/euphs$Yc_sd
euphs$Anomaly_yr <- (euphs$Abundance - meanNaN(euphs$Abundance)) / sdNaN(euphs$Abundance)

# fill region and taxa so they arent NA
euphs <- euphs %>%
  fill(Region, taxa, .direction = "downup") 

# give abundance data day & month just to match w index
euphs$date <- as.Date(paste0(euphs$Year, "-03-01"))

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

# z score PDO
PDO$pdo_z <- (PDO$pdo - meanNaN(PDO$pdo)) / sdNaN(PDO$pdo)

meanNaN(PDO$pdo_z)  # should be close to 0
sdNaN(PDO$pdo_z)    # should be close to 1
hist(PDO$pdo_z) 
plot(PDO$time, PDO$pdo_z, type = "l")

## ------------------------------------------ ##
#  Function for integrations (matching MATLAB/Manu's approach) -----
## ------------------------------------------ ##
# AR(1) style smoothing
recursive_integration <- function(x, tau, dt = 1) {
  # calculate autoregressive weight (how much past value is retained)
  # tau = memory in months, dt = time step 1 month
  alpha <- 1 - dt / tau # autoregressive weight
  y <- numeric(length(x))
  y[1] <- x[1] # initial condition
  for (i in 1:(length(x) - 1)) {
    y[i + 1] <- alpha * y[i] + x[i] * dt
  }
  return(y)
}

# apply integration
tau_bio <- 24        # biology memory in months
dt <- 1              # monthly time steps

PDO$pdo_int <- recursive_integration(PDO$pdo_z, tau = tau_bio, dt = dt)

# zscore after integration
PDO$pdo_int_z <- (PDO$pdo_int - meanNaN(PDO$pdo_int)) / sdNaN(PDO$pdo_int) #PDO$bio_z <- scale(PDO$bio)[,1]  
PDO$time <- as.Date(PDO$time)

## interpolate (bio is yearly, pdo is monthly)
euphs$pdo_orig <- approx(PDO$time, PDO$pdo_z, xout = euphs$date)$y
euphs$pdo_int <- approx(PDO$time, PDO$pdo_int_z, xout = euphs$date)$y

# correlations
cor_raw <- cor(euphs$Anomaly_yr, euphs$pdo_orig)
cor_int <- cor(euphs$Anomaly_yr, euphs$pdo_int)

# Get p-values
pval_raw <- cor.test(euphs$Anomaly_yr, euphs$pdo_orig)$p.value
pval_int <- cor.test(euphs$Anomaly_yr, euphs$pdo_int)$p.value


## export 
# save euphausiid data 
write.csv(euphs[, c("date", "Anomaly_yr")], "output/CCE/nsimplex_bio_output.csv", row.names = FALSE)

# save PDO monthly z-scored
write.csv(PDO[, c("time", "pdo_z")], "output/CCE/pdo_z_output.csv", row.names = FALSE)





#cor(euphs$Anomaly_yr, euphs$pdo_driver)
#cor.test(euphs$Anomaly_yr, euphs$pdo_driver)




ggplot(PDO, aes(x = time)) +
  geom_line(aes(y = pdo_z), color = "black", linewidth = 1.2) +
  geom_line(aes(y = bio_z), color = "red", linewidth = 1.2) +
  labs(y = "Standardized Value", x = "Year", title = "PDO (z-scored) and 2× Integrated Signal") +
  theme_minimal()

### to test sensitivity 
taus <- seq(6, 120, by = 6)  # test from 0.5 years to 10 years in 0.5-year steps 
#taus <- c(6, 12, 24)

cor_vals <- numeric(length(taus))

for (i in seq_along(taus)) {
  tau <- taus[i]  # current damping timescale in months
  
  # Apply the recursive integration using your defined function
  # This integrates the standardized PDO time series with memory = tau
  pdo_int <- recursive_integration(PDO$pdo_z, tau = tau, dt = 1)
  
  # Z-score the integrated PDO time series (so it's centered/scaled before correlation)
  pdo_int_z <- (pdo_int - mean(pdo_int, na.rm = TRUE)) / sd(pdo_int, na.rm = TRUE)
  
  # Interpolate the integrated PDO (monthly) to match the yearly euphausiid dates
  pdo_interp <- approx(PDO$time, pdo_int_z, xout = euphs$date)$y
  
  # Calculate correlation between yearly euphausiid anomalies and the integrated PDO
  cor_vals[i] <- cor(euphs$Anomaly_yr, pdo_interp, use = "pairwise.complete.obs")
}

# Plot the results: Correlation as a function of timescale (in years)
plot(taus / 12, cor_vals, type = "l", lwd = 2, col = "black", 
     xlab = "Timescale τ (years)", 
     ylab = "Correlation",
     main = "Sensitivity to τ (biological memory)")



for (tau in taus) {
  PDO[[paste0("bio_tau", tau)]] <- recursive_integration(PDO$ocean, tau = tau)
}


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

calculateIntegrations = function(data, tau = 365*2, 
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

pdoInt <- calculateIntegrations(PDO) #if your R window is too small may give err
pdoInt$time <- as.Date(pdoInt$time)
#write.csv(pdoInt, "output/CCE/PDO_CCE_integrated.csv")

## ------------------------------------------ ##
#     Correlations -----
## ------------------------------------------ ##
# 3) Test correlations between bio TS + each driver TS 

driver1 <- pdoInt[, c(1,3)]
int1 <- pdoInt[, c(1,4)]
#int2 <- pdoInt[, c(1,5)]

# interpolate driver values - matching dates/times 
interpolated_driver <- approx(driver1$time, driver1[[2]], xout = euphs$date)$y
interpolated_int1 <- approx(int1$time, int1[[2]], xout = euphs$date)$y
#interpolated_int2 <- approx(int2$time, int2[[2]], xout = euphs$date)$y

# calculate correlations
cor1 <- cor(euphs$Anomaly_yr, interpolated_driver)
cor2 <- cor(euphs$Anomaly_yr, interpolated_int1)
#cor3 <- cor(euphs$Anomaly_yr, interpolated_int2)

# get conventional pvalues 
pval1 <- cor.test(euphs$Anomaly_yr, interpolated_driver)$p.value
pval2 <- cor.test(euphs$Anomaly_yr, interpolated_int1)$p.value
#pval3 <- cor.test(euphs$Anomaly_yr, interpolated_int2)$p.value

wcor_CCE <- data.frame(site = "CCE",
                       taxa = "Nsimplex",
                       season = "spring",
                       cor_Norm = cor1, 
                       cor_Int = cor2, 
                       #cor_DouInt = cor_DouInt, 
                       pval_Norm = pval1, 
                       pval_Int = pval2
                       #pval_DouInt = pval_DouInt
)
wcor_CCE
#write.csv(wcor_CCE, "output/CCE/weightedCor_CCE.csv")

## ------------------------------------------ ##
#     Plots -----
## ------------------------------------------ ## 
(NormPlot <- ggplot() +
  # Driver time series (PDO normalized)
  geom_line(data = driver1, aes(x = time, y = pdoNorm, color = "PDO"), 
            linewidth = 1.2) +
  # Nsimplex time series (Abundance)
  geom_line(data = euphs, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
            linewidth = 1.2) +
  # annotate correlation coefficient
  annotate("text", x = as.Date("1960-01-01"), y = 3,
           label = sprintf("r == %.2f", cor1),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  # annotate p-value
  annotate("text", x = as.Date("1960-01-01"), y = 2.4,
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
  geom_line(data = int1, aes(x = time, y = pdoInt, color = "PDO"), 
            linewidth = 1.2) +
  # Nsimplex time series (Abundance)
  geom_line(data = euphs, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
            linewidth = 1.2) +
  # annotate correlation coefficient
  annotate("text", x = as.Date("1960-01-01"), y = 3,
           label = sprintf("r == %.2f", cor2),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  # annotate p-value
  annotate("text", x = as.Date("1960-01-01"), y = 2.4,
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
ggsave("figures/CCE/conventional_CCE_Ns_spring.png", 
       final_plot, 
       width = 8, height = 8, dpi = 300, 
       bg = "white")
##################################################################


Int2Plot <- ggplot() +
  # Driver time series (Integrated PDO)
  geom_line(data = int2, aes(x = time, y = pdoDInt, color = "MEI"), 
            linewidth = 1.2) +
  # Nsimplex time series (Abundance)
  geom_line(data = euphs, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
            linewidth = 1.2) +
  # annotate correlation coefficient
  annotate("text", x = as.Date("1960-01-01"), y = 3,
           label = sprintf("r == %.2f", cor3),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
  # annotate p-value
  annotate("text", x = as.Date("1960-01-01"), y = 2.4,
           label = if (pval3 < 0.001) {
             as.expression(bquote(italic(p) < 0.001))
           } else {
             as.expression(bquote(italic(p) == .(round(pval3, 3))))
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
  scale_x_date(breaks = seq(as.Date("1950-01-01"), as.Date("2022-01-01"), by = "10 years"),
               date_labels = "%Y", expand = c(0, 0)) +
  coord_cartesian(xlim = as.Date(c("1950-01-01", "2022-01-01"))) +
  # manual color legend
  scale_color_manual(values = c("PDO" = "red", "Abundance" = "blue")) +
  # y-axis label
  labs(y = "Double Integrated PDO")

Int2Plot <- Int2Plot + theme(legend.position = "none")