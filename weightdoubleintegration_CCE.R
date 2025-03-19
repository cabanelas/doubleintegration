################################################################################
#############          Pelagic Synthesis           #############################
#############             FEB-2025                 #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas 
################################################################################
##### CCE
##### species?
##### season?
##### 2 year integration time 
##### driver???

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

abundance_data <- read.csv(file.path("raw",
                                     "CCE",
                                     "Euphausiids_CCE.csv"))

## calculate anomalies
abundance_data$meanAbu <- mean(abundance_data$abundance)
abundance_data$stdAbu <- sd(abundance_data$abundance)
abundance_data$Anomaly_yr2 <- (abundance_data$abundance - abundance_data$meanAbu)/abundance_data$stdAbu
##

# year column to numeric
abundance_data$Year <- as.numeric(abundance_data$Year)

# Add 3.5/12 to each year (following M.Stukel Matlab script)
abundance_data$Year <- abundance_data$Year + 3.5/12
# OR abundance_data_a$date <- as.Date(paste0(abundance_data_a$Year, "-03-01"))

abundance_data <- abundance_data[, c(1,8)]
abundance_data_matrix <- as.matrix(abundance_data1)

abundance_data <- abundance_data %>%
  rename(date = Year)

## ------------------------------------------ ##
#            Driver Data -----
## ------------------------------------------ ##

pdo <- read.csv(file.path("raw",
                   "CCE",
                   "PDO.csv"))

pdo1 <-  pdo %>%
  pivot_longer(cols = Jan:Dec,
               names_to = "month",
               values_to = "pdo")

pdo1$Date <- zoo::as.yearmon(paste(pdo1$Year, pdo1$month), "%Y %b")

pdo1 <- pdo1 %>%
  rename(time = Date)

pdo1$time <- as.Date(paste(pdo1$time, "01"), 
                     format = "%b %Y %d")

pdo1 <- pdo1 %>% filter(
  #Year> 1944 & 
  pdo < 99)

pdo1 <- pdo1 %>% select(time, pdo)

pdo1 <- pdo1 %>% mutate(time = as.POSIXct(time, 
                                          format = "%Y-%m-%d")) 

pdo1 <- as.data.frame(pdo1)
pdo1$time <- pdo1$time - months(1)#check... shifting month back 1
#pdo1$date <- as.Date(paste0(pdo1$Year, "-03-01"))


## ------------------------------------------ ##
#     Function for integrations -----
## ------------------------------------------ ##
# 2) Driver TS 
#     a. Raw TS (at higher resolution than bio TS)
#     b. Temporally average (backwards in time) over characteristic 
#        time span of organism - 1st integration
#     c. Temporally average again - 2nd integration 

calculateIntegrations = function(data, tau = 365*2, f = function(x) {mean(x, na.rm = T)}) {
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

pdoInt1 <- calculateIntegrations(pdo1)
pdoInt1$time <- as.Date(pdoInt1$time)


## ------------------------------------------ ##
#     Correlations -----
## ------------------------------------------ ##
# 3) Test correlations between bio TS + each driver TS 

driver1 <- pdoInt1[, c(1,3)]
int1 <- pdoInt1[, c(1,4)]
int2 <- pdoInt1[, c(1,5)]

interpolated_driver <- approx(driver1$time, driver1[[2]], xout = abundance_data$date)$y
interpolated_int1 <- approx(int1$time, int1[[2]], xout = abundance_data$date)$y
interpolated_int2 <- approx(int2$time, int2[[2]], xout = abundance_data$date)$y

cor1 <- cor(abundance_data$Anomaly_yr, interpolated_driver)
cor2 <- cor(abundance_data$Anomaly_yr, interpolated_int1)
cor3 <- cor(abundance_data$Anomaly_yr, interpolated_int2)

pval1 <- cor.test(abundance_data$Anomaly_yr, interpolated_driver)$p.value
pval2 <- cor.test(abundance_data$Anomaly_yr, interpolated_int1)$p.value
pval3 <- cor.test(abundance_data$Anomaly_yr, interpolated_int2)$p.value









## ------------------------------------------ ##
#     Plots -----
## ------------------------------------------ ## 



driver1 <- driver1[driver1$time >= as.Date("1950-01-01"), ]
int1 <- int1[int1$time >= as.Date("1950-01-01"), ]
int2 <- int2[int2$time >= as.Date("1950-01-01"), ]

NormPlotb <- ggplot() +
  geom_line(data = driver1, aes(x = time, y = pdoNorm), 
            color = "red", size = 1.2) +
  geom_line(data = abundance_data, aes(x = date, y = Anomaly_yr), 
            color = "blue", size = 1.2) +
  annotate("text", x = as.Date("1953-01-01"), y = 3,
           label = sprintf("rho == %.4f", cor1),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 4.5) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks.x = element_blank()) + 
  scale_x_date(breaks = seq(as.Date("1950-01-01"), as.Date("2020-01-01"), by = "10 years"),
               date_labels = "%Y",
               expand = c(0, 0))

Int1Plotb <- ggplot() +
  geom_line(data = int1, aes(x = time, y = pdoInt), 
            color = "red", size = 1.2) +
  geom_line(data = abundance_data_a, aes(x = date, y = Anomaly_yr), 
            color = "blue", size = 1.2) +
  annotate("text", x = as.Date("1953-01-01"), y = 3,
           label = sprintf("rho == %.4f", cor2),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 4.5) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks.x = element_blank()) + 
  scale_x_date(breaks = seq(as.Date("1950-01-01"), as.Date("2020-01-01"), by = "10 years"),
               date_labels = "%Y",
               expand = c(0, 0))

Int2Plotb <- ggplot() +
  geom_line(data = int2, aes(x = time, y = pdoDInt),
            color = "red", size = 1.2) +
  geom_line(data = abundance_data_a, aes(x = date, y = Anomaly_yr), 
            color = "blue", size = 1.2) +
  annotate("text", x = as.Date("1953-01-01"), y = 3,
           label = sprintf("rho == %.4f", cor3),
           parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 4.5) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14, color = "black")) + 
  scale_x_date(breaks = seq(as.Date("1950-01-01"), as.Date("2020-01-01"), by = "10 years"),
               date_labels = "%Y",
               expand = c(0, 0))

# Combine the three plots into one grid
plot_grid(NormPlotb, Int1Plotb, Int2Plotb, 
          ncol = 1, align = "v")