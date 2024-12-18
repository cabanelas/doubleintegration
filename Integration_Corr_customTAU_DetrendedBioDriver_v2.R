################################################################################
#############          Pelagic Synthesis           #############################
#############             MAR-2024                 #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas 
################################################################################
## "Single" Integration Analysis NES - ECOMON - since im using oceanic variables
# Script #4 : Integration_Corr_customTAU_DetrendedBioDriver

# script to perform integrations, correlations, and plot
# code allows to specify custom TAU (generation/damping time scale) for diff taxa 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##

library(tidyverse) #v2.0.0
library(gridExtra) #v2.3 
library(cowplot) #v1.1.3

#install.packages("devtools")
#remotes::install_github("noaa-edab/ecodata",build_vignettes=TRUE) 
#or pak::pkg_install("noaa-edab/ecodata")
#library(ecodata) #v5.0.1
#https://github.com/NOAA-EDAB/ecodata

library(listviewer) #v4.0.0; for looking at lists interactively 
#library(magrittr) #v2.0.3; map_dfr
library(here) #v1.0.1; easily build path to files  
################################################################################

## ------------------------------------------ ##
#            Biology Data -----
## ------------------------------------------ ##

# can download EcoMon data from 
#https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0187513 
# NCEI Accession 0187513 v3.3 

# CREATED IN calculateARcoefficient_bio_v3
zscore <- read.csv(file.path("output",
                             "detrended_BIOLOGY_time_series.csv"))

ggplot(zscore, aes(x=year, y=Anomaly_yr)) + 
  geom_line() + 
  facet_grid(taxa ~ season+Region)

# check SD and mean of anomalies by group
zscore %>%
  group_by(taxa, Region, season) %>%
  summarize(mean_Anomaly_yr = mean(Anomaly_yr, na.rm = TRUE),
            sd_Anomaly_yr = sd(Anomaly_yr, na.rm = TRUE)) %>%
  print(n = 50)

# remove ammspp fall SNE since == 0 
zscore <- zscore[complete.cases(zscore$Anomaly_yr), ]

###############################################################################
# 2) Driver TS 
#     a. Raw TS (at higher resolution than bio TS)
#     b. Temporally average (backwards in time) over characteristic 
#        time span of organism - 1st integration
#     c. Temporally average again - 2nd integration 

#standardize = function(x) {
#  (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
#}

## ------------------------------------------ ##
#     Function for integrations -----
## ------------------------------------------ ##
calculateIntegrations = function(data, tau, f = function(x) {mean(x, na.rm = T)}) {
  for (n in names(data)[-1]) {
    data[[paste0(n,'Norm')]] = (data[[n]] - mean(data[[n]], na.rm = T)) / sd(data[[n]]) #norm
    
    data[[paste0(n,'Int')]] = NA
    
    for (i in 2:nrow(data)) {
      k = data[,1] <= data[i,1] & data[,1] > data[i,1] - tau * 86400 #secs
      data[[paste0(n,'Int')]][i] = f(data[[paste0(n,'Norm')]][k])
    }
  }
  data
}

# tau values for each taxa
tau_values <- c(
  "ctyp" = 60, #2 months
  "calfin" = 240, #8 months
  "pseudo" = 90 #3 months
)


#There are 86,400 seconds in a day (60 seconds/minute * 60 minutes/hour * 24 hours/day).
#multiplying tau by 86,400 gives the equivalent time window in seconds.
################################################################################
################################################################################

## Drivers - AMO, NAO, AO, GSI
#                               Driver Data 
##### ALL THESE WERE DETRENDED IN: calculateARcoefficient_drivers_v3 script
#by fitting linear regression
################################################################################

##                  AMO - Atlantic Multidecadal Oscillation                   ##
# original came from: https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v5/index/ersst.v5.amo.dat
# 1973-2023

##                            NAO - North Atlantic Oscillation                ##
# original came from: https://www.ncei.noaa.gov/access/monitoring/nao/
#https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii
# 1973-2023

#also available through this API download_nao() rsoi package
#data comes from "https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii.table"

##                           AO - Arctic Oscillation                          ##
# original came from: https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/monthly.ao.index.b50.current.ascii
# https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/ao.shtml
# 1973-2023

#also available through API download_ao() from package rsoi
#from "https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/monthly.ao.index.b50.current.ascii.table"

##                               GSI - Gulf Stream Index                      ##
# data(package="ecodata"); 
#gsi <- ecodata::gsi
#write.csv(gsi, "raw/gsi.csv") downloaded 26-MAR-2024
# 1954-2023
###################
## ------------------------------------------ ##
#            Driver Data -----
## ------------------------------------------ ##

drivers <- read.csv(file.path("output",
                              "detrended_driver_time_series.csv"))

amo <- drivers %>% select(time, AMO_detrended) %>%
  rename_with(~ tolower(.), ends_with("_detrended")) %>%
  mutate(time = as.POSIXct(time, format = "%Y-%m-%d")) 
nao <- drivers %>% select(time, NAO_detrended) %>%
  rename_with(~ tolower(.), ends_with("_detrended")) %>%
  mutate(time = as.POSIXct(time, format = "%Y-%m-%d")) 
ao <- drivers %>% select(time, AO_detrended) %>%
  rename_with(~ tolower(.), ends_with("_detrended")) %>%
  mutate(time = as.POSIXct(time, format = "%Y-%m-%d")) 
gsi <- drivers %>% select(time, GSI_detrended) %>%
  rename_with(~ tolower(.), ends_with("_detrended")) %>%
  mutate(time = as.POSIXct(time, format = "%Y-%m-%d")) 

# Create a list of data frames, one for each tau value
amoInt_list <- lapply(tau_values, function(tau) calculateIntegrations(amo, tau))
naoInt_list <- lapply(tau_values, function(tau) calculateIntegrations(nao, tau))
aoInt_list <- lapply(tau_values, function(tau) calculateIntegrations(ao, tau))
gsiInt_list <- lapply(tau_values, function(tau) calculateIntegrations(gsi, tau))

listviewer::jsonedit(amoInt_list)
listviewer::jsonedit(naoInt_list)

###############
# plot TS of drivers with 3 double integration 

plotData <- function(data, int_list, tau_values) {
  par(mfrow = c(2,1))
  for (n in names(data)[-1]) {
    for (tau_name in names(int_list)) {
      driver_name <- names(int_list[[tau_name]])[2] #Extract driver name for title
      plot(data[,1], int_list[[tau_name]][[n]], 
           type = 'l', ylab = n, xlab = 'Date') #Plot each integration
      title(paste(driver_name, "  TAU:", tau_values[tau_name]), 
            outer = TRUE, line = -1) # Add driver name and TAU value to title
      plot(data[,1], int_list[[tau_name]][[paste0(n,'Int')]], 
           type = 'l', ylab = paste0(n, ' Int'), xlab = 'Date') # Plot Int
    }
  }
}

# plot each driver with varying tau values 
plotData(amo, amoInt_list, tau_values)
plotData(nao, naoInt_list, tau_values)
plotData(ao, aoInt_list, tau_values)
plotData(gsi, gsiInt_list, tau_values)
par(mfrow = c(1,1))
################################################################################

## ------------------------------------------ ##
#     Correlations -----
## ------------------------------------------ ##

# Test correlations between bio TS + each driver TS 

# List to store the correlations and pvalues for each combination/group
correlation_list <- list()
p_value_list <- list()

# Get unique combinations of season, taxa, and Region in zscore
combinations <- unique(zscore[, c("season", "taxa", "Region")])

# Iterate over each unique combination
for (i in seq_len(nrow(combinations))) {
  # Extract the current combination
  season <- combinations[i, "season"]
  taxa <- combinations[i, "taxa"]
  Region <- combinations[i, "Region"]
  
  cat("Processing group:", paste(season, taxa, Region, sep = "_"), "\n")
  
  # Filter data for the current combination
  group_data <- zscore[zscore$season == season & 
                       zscore$taxa == taxa &
                       zscore$Region == Region, ]
  
  # Convert date to Date format
  group_data$date <- as.Date(group_data$date)
  
  # List to store data for each combination of taxa-group and driver
  group_correlations <- list()
  group_pvalues <- list()
  
  # Iterate over each driver
  for (driver_name in c("amo", "nao", "ao", "gsi")) {
    cat("Calculating correlations for driver:", driver_name, "\n")
    
    # Extract driver data 
    driver_data_list <- switch(driver_name,
                               "amo" = amoInt_list,
                               "nao" = naoInt_list,
                               "ao" = aoInt_list,
                               "gsi" = gsiInt_list)
    
    # Iterate over each taxa and its corresponding tau value (assuming tau_values is defined somewhere)
    for (taxa_name in names(tau_values)) {
      # Check if the taxa_name matches the current combination
      if (taxa == taxa_name) {
        # Get the driver data for the current taxa
        driver_data <- driver_data_list[[taxa_name]]
        
        # Convert driver_data$time to Date format
        driver_data$time <- as.Date(driver_data$time)
        
        # Interpolate driver values
        interpolated_driver <- approx(driver_data$time, driver_data[[3]], xout = group_data$date)$y
        interpolated_int1 <- approx(driver_data$time, driver_data[[4]], xout = group_data$date)$y
        
        # Calculate correlations
        correlations <- c(
          cor(group_data$Anomaly_yr, interpolated_driver, use = "complete.obs"),
          cor(group_data$Anomaly_yr, interpolated_int1, use = "complete.obs")
        )
        
        p_values <- c(
          cor.test(group_data$Anomaly_yr, interpolated_driver)$p.value,
          cor.test(group_data$Anomaly_yr, interpolated_int1)$p.value
        )

        # Store correlations, p-values, and Bonferroni-adjusted p-values
        group_correlations[[paste0(season, "_", taxa, "_", Region, "_", driver_name, "_", taxa_name)]] <- correlations
        group_pvalues[[paste0(season, "_", taxa, "_", Region, "_", driver_name, "_", taxa_name)]] <- p_values
      }
    }
  }
  
  # Store correlations and p-values for the current group combination
  correlation_list[[paste(season, taxa, Region, sep = "_")]] <- group_correlations
  p_value_list[[paste(season, taxa, Region, sep = "_")]] <- group_pvalues
}

correlation_list
p_value_list

map(correlation_list, as.data.frame) %>%
  map(t)

lapply(correlation_list, as.data.frame)

# create df with correlation- and p-values
#empty vectors to store data
names <- c()
cor_Norm <- c()
cor_Int <- c()
#cor_DouInt <- c()
pval_Norm <- c()
pval_Int <- c()
#pval_DouInt <- c()

# Loop through the list to extract names, correlation values, and p-values
for (season in names(correlation_list)) {
  for (taxa in names(correlation_list[[season]])) {
    names <- c(names, taxa)
    cor_Norm <- c(cor_Norm, correlation_list[[season]][[taxa]][1])
    cor_Int <- c(cor_Int, correlation_list[[season]][[taxa]][2])
    #cor_DouInt <- c(cor_DouInt, correlation_list[[season]][[taxa]][3])
    pval_Norm <- c(pval_Norm, p_value_list[[season]][[taxa]][1])
    pval_Int <- c(pval_Int, p_value_list[[season]][[taxa]][2])
    #pval_DouInt <- c(pval_DouInt, p_value_list[[season]][[taxa]][3])
  }
}
# Create df
df <- data.frame(name = names, 
                 cor_Norm = cor_Norm, 
                 cor_Int = cor_Int, 
                 #cor_DouInt = cor_DouInt, 
                 pval_Norm = pval_Norm, 
                 pval_Int = pval_Int
                 #pval_DouInt = pval_DouInt
                )

df <- df %>%
  separate(name, 
           into = c("taxa", "region", "season", "driver", "ctyp"), 
           sep = "_", 
           remove = FALSE)

df <- df %>%
  select(-c(name, ctyp))

#write.csv(df, "output/IntCorrelations_DetrendBioDriver_6090240TAU.csv")


## ------------------------------------------ ##
#     PLOTS -----
## ------------------------------------------ ##

# List to store the final plots 
final_plots_list <- list()

output_dir <- here("figures/21JUL_DetrendDriverBIO_6090240tau")

driver_full_names <- list(
  amo = "Atlantic Multidecadal Oscillation",
  nao = "North Atlantic Oscillation",
  ao = "Arctic Oscillation",
  gsi = "Gulf Stream Index"
)

# Define full taxa names
full_taxa_names <- list(
  ctyp = "Centropages typicus",
  calfin = "Calanus finmarchicus",
  pseudo = "Pseudocalanus spp."
)

# Define full region names
full_region_names <- list(
  "GOM" = "Gulf of Maine",
  "MAB" = "Mid-Atlantic Bight",
  "GB" = "Georges Bank",
  "SNE" = "Southern New England"
)

# Iterate over each group in correlation_list
for (group_name in names(correlation_list)) {
  cat("Processing group:", group_name, "\n")
  
  # Extract abundance data for the current group
  abundance_data <- zscore[zscore$season == strsplit(group_name, "_")[[1]][1] &
                             zscore$taxa == strsplit(group_name, "_")[[1]][2] &
                             zscore$Region == strsplit(group_name, "_")[[1]][3], ]
  
  # List to store the final plots for each driver
  group_plots_list <- list()

  # Iterate over each atmospheric driver
  for (driver_name in c("amo", "nao", "ao", "gsi")) {  # Drivers
  cat("Calculating correlations for driver:", driver_name, "\n")
  
  # Extract correlations for the current driver
  correlations <- correlation_list[[group_name]]
  pvalues <- p_value_list[[group_name]]
  
  # Extract the correlation values for the current driver, taxa, and group combination
  current_correlations <- correlations[[paste0(group_name, "_", 
                                               driver_name, "_", 
                                               strsplit(group_name, "_")[[1]][2])]]
  current_pvals <- pvalues[[paste0(group_name, "_", 
                                   driver_name, "_", 
                                   strsplit(group_name, "_")[[1]][2])]]
  
  # Extract driver data
  driver_data_list <- switch(driver_name,
                             "amo" = amoInt_list,
                             "nao" = naoInt_list,
                             "ao" = aoInt_list,
                             "gsi" = gsiInt_list)
  
  # Get driver data for the current taxa
  driver_data <- driver_data_list[[strsplit(group_name, "_")[[1]][2]]]
  
  # Extract taxa, Region, and season for the current group
  taxa <- abundance_data$taxa  
  Region <- abundance_data$Region  
  season <- abundance_data$season
  
  taxa1 <- as.character(full_taxa_names[taxa])
  region1 <- as.character(full_region_names[Region])
  
  driver_label <- toupper(driver_name)
  driver_legend <- paste(driver_full_names[[driver_name]])
  
  abundance_data$date <- as.POSIXct(abundance_data$date)
  driver_data$time <- as.POSIXct(driver_data$time)
  
  # plots 
  NormPlot <- ggplot() +
    geom_line(data = driver_data, aes(x = time, 
                                      y = .data[[paste0(driver_name, "_detrendedNorm")]],
                                      color = driver_legend), 
               size = 1.2) +
    geom_smooth(data = driver_data, aes(x = time, 
                                        y = .data[[paste0(driver_name, "_detrendedNorm")]]), 
                method = "loess", 
                span = 0.1, 
                se = FALSE, 
                size = 1,
                color = "black") +
    geom_line(data = abundance_data, aes(x = date, 
                                         y = Anomaly_yr,
                                         color = "Abundance"), 
              size = 1.2) +
    annotate("text", x = as.POSIXct("1975-01-01"), y = 3,
             label = paste("r =", format(round(current_correlations[1], 2), 
                                         nsmall = 2)),
             parse = F, hjust = 0, vjust = 1, color = "black", size = 5.5) +
    annotate("text", x = as.POSIXct("1975-01-01"), y = 2.4,
             label = if (current_pvals[1] < 0.001) {
               as.expression(bquote(italic(P) < 0.001))
             } else {
               as.expression(bquote(italic(P) == .(round(current_pvals[1], 3))))
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
    scale_x_datetime(breaks = seq(as.POSIXct("1972-01-01"), 
                                  as.POSIXct("2024-01-01"), 
                                  by = "10 years"),
                     date_labels = "%Y", expand = c(0, 0)) +
    coord_cartesian(xlim = as.POSIXct(c("1972-01-01", "2024-01-01"))) +
    scale_color_manual(values = setNames(c("red", "blue"), 
                                         c(driver_legend, "Abundance"))) +
    labs(y = driver_label)
  
  Int1Plot <- ggplot() +
    geom_line(data = driver_data, aes(x = time, 
                                      y = .data[[paste0(driver_name, "_detrendedInt")]], 
                                      color = driver_legend), 
              size = 1.2) +
    geom_line(data = abundance_data, aes(x = date, y = Anomaly_yr, 
                                         color = "Abundance"), size = 1.2) +
    annotate("text", x = as.POSIXct("1975-01-01"), y = 3,
             label = paste("r =", format(round(current_correlations[2], 2), 
                                         nsmall = 2)),
             parse = F, hjust = 0, vjust = 1, color = "black", size = 5.5) +
    annotate("text", x = as.POSIXct("1975-01-01"), y = 2.4,
             label = if (current_pvals[2] < 0.001) {
               as.expression(bquote(italic(P) < 0.001))
             } else {
               as.expression(bquote(italic(P) == .(round(current_pvals[2], 3))))
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
    scale_x_datetime(breaks = seq(as.POSIXct("1972-01-01"), 
                                  as.POSIXct("2024-01-01"), 
                                  by = "10 years"),
                     date_labels = "%Y", expand = c(0, 0)) +
    coord_cartesian(xlim = as.POSIXct(c("1972-01-01", "2024-01-01"))) +
    scale_color_manual(values = setNames(c("red", "blue"), 
                                         c(driver_legend, "Abundance"))) +
    labs(y = paste("Integrated", driver_label))
  
  legend <- cowplot::get_plot_component(NormPlot, 'guide-box-bottom', 
                                        return_all = TRUE)
  cowplot::ggdraw(legend)
  
  # Remove legend
  NormPlot <- NormPlot + theme(legend.position = "none")
  Int1Plot <- Int1Plot + theme(legend.position = "none")
  
  # Combine plots into one grid
  combined_plots <- plot_grid(NormPlot, Int1Plot, 
                              ncol = 1, align = "v")
  
  # Add plot title 
  title <- ggdraw() +
    draw_label(label = bquote(paste(italic(.(taxa1)), " in ",
                                    .(region1), " in ",
                                    .(season))), 
               size = 18) +
    theme(plot.background = element_rect(fill = "transparent"))
  
  # Combine title and plots
  final_plot <- plot_grid(title, combined_plots, legend, ncol = 1, 
                          rel_heights = c(0.07, 0.8, 0.07)) 
  
  file_name <- paste0(driver_name, "_", group_name, ".png")
  
  # Save plot
  ggsave(file.path(output_dir, file_name), 
         final_plot, 
         width = 8, height = 8, dpi = 300, 
         bg = "white")
  
  #right now only problem is that it saves it with border 
  #could try panel_border() or 
  #plot.background = element_rect(color = "white") inside theme()
  #but pretty sure the issue is with ggsave
  
  # Store the final plot for the current driver
  group_plots_list[[driver_name]] <- final_plot
  }
  # Store the final plots for the current group
  final_plots_list[[group_name]] <- group_plots_list
}

final_plots_list[["ctyp_MAB_summer"]]
final_plots_list[["calfin_SNE_fall"]]


## ------------------------------------------ ##
#     PLOTS ----- NO LOESS LINE
## ------------------------------------------ ##

# List to store the final plots 
final_plots_list <- list()

output_dir <- here("figures/21JUL_DetrendDriverBIO_6090240tau_noLOESS")

# Iterate over each group in correlation_list
for (group_name in names(correlation_list)) {
  cat("Processing group:", group_name, "\n")
  
  # Extract abundance data for the current group
  abundance_data <- zscore[zscore$season == strsplit(group_name, "_")[[1]][1] &
                             zscore$taxa == strsplit(group_name, "_")[[1]][2] &
                             zscore$Region == strsplit(group_name, "_")[[1]][3], ]
  
  # List to store the final plots for each driver
  group_plots_list <- list()
  
  # Iterate over each atmospheric driver
  for (driver_name in c("amo", "nao", "ao", "gsi")) {  # Drivers
    cat("Calculating correlations for driver:", driver_name, "\n")
    
    # Extract correlations for the current driver
    correlations <- correlation_list[[group_name]]
    pvalues <- p_value_list[[group_name]]
    
    # Extract the correlation values for the current driver, taxa, and group combination
    current_correlations <- correlations[[paste0(group_name, "_", 
                                                 driver_name, "_", 
                                                 strsplit(group_name, "_")[[1]][2])]]
    current_pvals <- pvalues[[paste0(group_name, "_", 
                                     driver_name, "_", 
                                     strsplit(group_name, "_")[[1]][2])]]
    
    # Extract driver data
    driver_data_list <- switch(driver_name,
                               "amo" = amoInt_list,
                               "nao" = naoInt_list,
                               "ao" = aoInt_list,
                               "gsi" = gsiInt_list)
    
    # Get driver data for the current taxa
    driver_data <- driver_data_list[[strsplit(group_name, "_")[[1]][2]]]
    
    # Extract taxa, Region, and season for the current group
    taxa <- abundance_data$taxa  
    Region <- abundance_data$Region  
    season <- abundance_data$season
    
    taxa1 <- as.character(full_taxa_names[taxa])
    region1 <- as.character(full_region_names[Region])
    
    driver_label <- toupper(driver_name)
    driver_legend <- paste(driver_full_names[[driver_name]])
    
    abundance_data$date <- as.POSIXct(abundance_data$date)
    driver_data$time <- as.POSIXct(driver_data$time)
    
    # plots 
    NormPlot <- ggplot() +
      geom_line(data = driver_data, aes(x = time, 
                                        y = .data[[paste0(driver_name, "_detrendedNorm")]],
                                        color = driver_legend), 
                size = 1.2) +
      geom_line(data = abundance_data, aes(x = date, 
                                           y = Anomaly_yr,
                                           color = "Abundance"), 
                size = 1.2) +
      annotate("text", x = as.POSIXct("1975-01-01"), y = 3,
               label = paste("r =", format(round(current_correlations[1], 2), 
                                           nsmall = 2)),
               parse = F, hjust = 0, vjust = 1, color = "black", size = 5.5) +
      annotate("text", x = as.POSIXct("1975-01-01"), y = 2.4,
               label = if (current_pvals[1] < 0.001) {
                 as.expression(bquote(italic(P) < 0.001))
               } else {
                 as.expression(bquote(italic(P) == .(round(current_pvals[1], 3))))
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
      scale_x_datetime(breaks = seq(as.POSIXct("1972-01-01"), 
                                    as.POSIXct("2024-01-01"), 
                                    by = "10 years"),
                       date_labels = "%Y", expand = c(0, 0)) +
      coord_cartesian(xlim = as.POSIXct(c("1972-01-01", "2024-01-01"))) +
      scale_color_manual(values = setNames(c("red", "blue"), 
                                           c(driver_legend, "Abundance"))) +
      labs(y = driver_label)
    
    Int1Plot <- ggplot() +
      geom_line(data = driver_data, aes(x = time, 
                                        y = .data[[paste0(driver_name, "_detrendedInt")]], 
                                        color = driver_legend), 
                size = 1.2) +
      geom_line(data = abundance_data, aes(x = date, y = Anomaly_yr, 
                                           color = "Abundance"), size = 1.2) +
      annotate("text", x = as.POSIXct("1975-01-01"), y = 3,
               label = paste("r =", format(round(current_correlations[2], 2), 
                                           nsmall = 2)),
               parse = F, hjust = 0, vjust = 1, color = "black", size = 5.5) +
      annotate("text", x = as.POSIXct("1975-01-01"), y = 2.4,
               label = if (current_pvals[2] < 0.001) {
                 as.expression(bquote(italic(P) < 0.001))
               } else {
                 as.expression(bquote(italic(P) == .(round(current_pvals[2], 3))))
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
      scale_x_datetime(breaks = seq(as.POSIXct("1972-01-01"), 
                                    as.POSIXct("2024-01-01"), 
                                    by = "10 years"),
                       date_labels = "%Y", expand = c(0, 0)) +
      coord_cartesian(xlim = as.POSIXct(c("1972-01-01", "2024-01-01"))) +
      scale_color_manual(values = setNames(c("red", "blue"), 
                                           c(driver_legend, "Abundance"))) +
      labs(y = paste("Integrated", driver_label))
    
    legend <- cowplot::get_plot_component(NormPlot, 'guide-box-bottom', 
                                          return_all = TRUE)
    cowplot::ggdraw(legend)
    
    # Remove legend
    NormPlot <- NormPlot + theme(legend.position = "none")
    Int1Plot <- Int1Plot + theme(legend.position = "none")
    
    # Combine plots into one grid
    combined_plots <- plot_grid(NormPlot, Int1Plot, 
                                ncol = 1, align = "v")
    
    # Add plot title 
    title <- ggdraw() +
      draw_label(label = bquote(paste(italic(.(taxa1)), " in ",
                                      .(region1), " in ",
                                      .(season))), 
                 size = 18) +
      theme(plot.background = element_rect(fill = "transparent"))
    
    # Combine title and plots
    final_plot <- plot_grid(title, combined_plots, legend, ncol = 1, 
                            rel_heights = c(0.07, 0.8, 0.07)) 
    
    file_name <- paste0(driver_name, "_", group_name, ".png")
    
    # Save plot
    ggsave(file.path(output_dir, file_name), 
           final_plot, 
           width = 8, height = 8, dpi = 300, 
           bg = "white")
    
    #right now only problem is that it saves it with border 
    #could try panel_border() or 
    #plot.background = element_rect(color = "white") inside theme()
    #but pretty sure the issue is with ggsave
    
    # Store the final plot for the current driver
    group_plots_list[[driver_name]] <- final_plot
  }
  # Store the final plots for the current group
  final_plots_list[[group_name]] <- group_plots_list
}
