################################################################################
#############          Pelagic Synthesis           #############################
#############                July 2024             #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas 
################################################################################

## Double Integration Analysis NES
# Script #6 : Final_plots_Conv_Boot_pvalues
# plotting everything together with both the conventional and bootstrapped pvals

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(tidyverse) #v2.0.0
library(here) #v1.0.1
library(cowplot) #v1.1.3

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##

# CREATED IN Boot_parallel_custom_TAU
res <- read.csv(file.path(
  #"output",
                             "results_boot_6090240TAU.csv")) %>%
  select(-X)


# CREATED IN calculateARcoefficient_drivers_v3
drivers <- read.csv(file.path("output",
                              "detrended_driver_time_series.csv"))


# CREATED IN calculateARcoefficient_bio_v3
zscore <- read.csv(file.path("output",
                             "detrended_BIOLOGY_time_series.csv")) %>%
  select(-X) %>%
  mutate(date = as.POSIXct(date, format = "%Y-%m-%d"))

zscore <- zscore[complete.cases(zscore$Anomaly_yr), ]


## ------------------------------------------ ##
#            Get driver TS for plotting -----
## ------------------------------------------ ##

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

# Create a list of data frames, one for each tau value
amoInt_list <- lapply(tau_values, function(tau) calculateIntegrations(amo, tau))
naoInt_list <- lapply(tau_values, function(tau) calculateIntegrations(nao, tau))
aoInt_list <- lapply(tau_values, function(tau) calculateIntegrations(ao, tau))
gsiInt_list <- lapply(tau_values, function(tau) calculateIntegrations(gsi, tau))

## ------------------------------------------ ##
#           Turn values into lists -----
## ------------------------------------------ ##

# Initialize correlation_list with the correct structure
correlation_list <- list()

# Loop through each row of the `res` data frame
for (i in 1:nrow(res)) {
  row <- res[i, ]
  
  # Extract relevant parts from the row for constructing the key
  group_name <- gsub("\\.", "_", row$group_name)  # Replace dots with underscores
  driver_name <- row$driver_name
  taxa_name <- row$taxa
  
  # Construct the key for the correlation_list
  key <- paste0(group_name, "_", driver_name, "_", taxa_name)
  
  # Replace any dots with underscores in the key
  key <- gsub("\\.", "_", key)
  
  # Check if the group already exists in the list
  if (!group_name %in% names(correlation_list)) {
    # Initialize an empty list for this group if it doesn't exist
    correlation_list[[group_name]] <- list()
  }
  
  # Check if the key already exists for this group
  if (!key %in% names(correlation_list[[group_name]])) {
    # Initialize an empty vector for this key if it doesn't exist
    correlation_list[[group_name]][[key]] <- c()
  }
  
  # Append the correlation values to the list entry
  correlation_list[[group_name]][[key]] <- c(correlation_list[[group_name]][[key]], 
                                             row$original_cor_norm, 
                                             row$original_cor_int)
}

# Print the result
print(correlation_list)



# Initialize p_value_list with the correct structure
p_value_list <- list()

# Loop through each row of the `res` data frame
for (i in 1:nrow(res)) {
  row <- res[i, ]
  
  # Extract relevant parts from the row for constructing the key
  group_name <- gsub("\\.", "_", row$group_name)  # Replace dots with underscores
  driver_name <- row$driver_name
  taxa_name <- row$taxa
  
  # Construct the key for the correlation_list
  key <- paste0(group_name, "_", driver_name, "_", taxa_name)

  # Replace any dots with underscores in the key
  key <- gsub("\\.", "_", key)
  
  # Check if the group already exists in the list
  if (!group_name %in% names(p_value_list)) {
    # Initialize an empty list for this group if it doesn't exist
    p_value_list[[group_name]] <- list()
  }
  
  # Check if the key already exists for this group
  if (!key %in% names(p_value_list[[group_name]])) {
    # Initialize an empty vector for this key if it doesn't exist
    p_value_list[[group_name]][[key]] <- c()
  }
  
  # Append the p-value values to the list entry
  p_value_list[[group_name]][[key]] <- c(p_value_list[[group_name]][[key]], row$original_pval_norm, row$original_pval_int)
}

# Print the result
print(p_value_list)



# Create a named list to store bootstrap p-values
bootstrap_p_val_list <- list()

# Loop through each row of the `res` data frame
for (i in 1:nrow(res)) {
  row <- res[i, ]
  
  # Extract relevant parts from the row for constructing the key
  group_name <- gsub("\\.", "_", row$group_name)  # Replace dots with underscores
  driver_name <- row$driver_name
  taxa_name <- row$taxa
  
  # Construct the key for the bootstrap_p_val_list
  key <- paste0(group_name, "_", driver_name, "_", taxa_name)
  
  # Replace any dots with underscores in the key
  key <- gsub("\\.", "_", key)
  
  # Check if the group already exists in the list
  if (!group_name %in% names(bootstrap_p_val_list)) {
    # Initialize an empty list for this group if it doesn't exist
    bootstrap_p_val_list[[group_name]] <- list()
  }
  
  # Check if the key already exists for this group
  if (!key %in% names(bootstrap_p_val_list[[group_name]])) {
    # Initialize an empty vector for this key if it doesn't exist
    bootstrap_p_val_list[[group_name]][[key]] <- c()
  }
  
  # Append the bootstrap p-value values to the list entry
  bootstrap_p_val_list[[group_name]][[key]] <- c(bootstrap_p_val_list[[group_name]][[key]], 
                                                 row$pval_pair, 
                                                 row$pval_pair_integrated)
}

# Print the result
print(bootstrap_p_val_list)

## ------------------------------------------ ##
#     PLOTS -----
## ------------------------------------------ ##

# List to store the final plots 
final_plots_list <- list()


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


output_dir <- here("figures/FINAL_boot_conv_pvals")
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
    boot_pvalues <- bootstrap_p_val_list[[group_name]]
    
    # Extract the correlation values for the current driver, taxa, and group combination
    current_correlations <- correlations[[paste0(group_name, "_", 
                                                 driver_name, "_", 
                                                 strsplit(group_name, "_")[[1]][2])]]
    current_pvals <- pvalues[[paste0(group_name, "_", 
                                     driver_name, "_", 
                                     strsplit(group_name, "_")[[1]][2])]]
    current_boot_pvals <- boot_pvalues[[paste0(group_name, "_", 
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
      annotate("text", x = as.POSIXct("1975-01-01"), y = 3.2,
               label = paste("r =", format(round(current_correlations[1], 2), 
                                           nsmall = 2)),
               parse = F, hjust = 0, vjust = 1, color = "black", size = 5.5) +
      annotate("text", x = as.POSIXct("1975-01-01"), y = 2.8,
               label = if (current_pvals[1] < 0.001) {
                 as.expression(bquote(italic(P)["conventional"] < 0.001))
               } else {
                 as.expression(bquote(italic(P)["conventional"] == .(round(current_pvals[1], 3))))
               },
               parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
      annotate("text", x = as.POSIXct("1975-01-01"), y = 2.4,
               label = if (current_boot_pvals[1] < 0.001) {
                 as.expression(bquote(italic(P)[bootstrap] < 0.001))
               } else {
                 as.expression(bquote(italic(P)[bootstrap] == .(round(current_boot_pvals[1], 3))))
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
      annotate("text", x = as.POSIXct("1975-01-01"), y = 3.2,
               label = paste("r =", format(round(current_correlations[2], 2), 
                                           nsmall = 2)),
               parse = F, hjust = 0, vjust = 1, color = "black", size = 5.5) +
      annotate("text", x = as.POSIXct("1975-01-01"), y = 2.8,
               label = if (current_pvals[2] < 0.001) {
                 as.expression(bquote(italic(P)["conventional"] < 0.001))
               } else {
                 as.expression(bquote(italic(P)["conventional"] == .(round(current_pvals[2], 3))))
               },
               parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
      annotate("text", x = as.POSIXct("1975-01-01"), y = 2.4,
               label = if (current_boot_pvals[2] < 0.001) {
                 as.expression(bquote(italic(P)[bootstrap] < 0.001))
               } else {
                 as.expression(bquote(italic(P)[bootstrap] == .(round(current_boot_pvals[2], 3))))
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

    # Store the final plot for the current driver
    group_plots_list[[driver_name]] <- final_plot
  }
  # Store the final plots for the current group
  final_plots_list[[group_name]] <- group_plots_list
}


## ------------------------------------------ ##
#     PLOTS with LOESS -----
## ------------------------------------------ ##
# List to store the final plots 
final_plots_list <- list()

output_dir <- here("figures/FINAL_boot_conv_pvals_withLOESS")

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
    boot_pvalues <- bootstrap_p_val_list[[group_name]]
    
    # Extract the correlation values for the current driver, taxa, and group combination
    current_correlations <- correlations[[paste0(group_name, "_", 
                                                 driver_name, "_", 
                                                 strsplit(group_name, "_")[[1]][2])]]
    current_pvals <- pvalues[[paste0(group_name, "_", 
                                     driver_name, "_", 
                                     strsplit(group_name, "_")[[1]][2])]]
    current_boot_pvals <- boot_pvalues[[paste0(group_name, "_", 
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
      annotate("text", x = as.POSIXct("1975-01-01"), y = 3.2,
               label = paste("r =", format(round(current_correlations[1], 2), 
                                           nsmall = 2)),
               parse = F, hjust = 0, vjust = 1, color = "black", size = 5.5) +
      annotate("text", x = as.POSIXct("1975-01-01"), y = 2.8,
               label = if (current_pvals[1] < 0.001) {
                 as.expression(bquote(italic(P)["conventional"] < 0.001))
               } else {
                 as.expression(bquote(italic(P)["conventional"] == .(round(current_pvals[1], 3))))
               },
               parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
      annotate("text", x = as.POSIXct("1975-01-01"), y = 2.4,
               label = if (current_boot_pvals[1] < 0.001) {
                 as.expression(bquote(italic(P)[bootstrap] < 0.001))
               } else {
                 as.expression(bquote(italic(P)[bootstrap] == .(round(current_boot_pvals[1], 3))))
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
      geom_smooth(data = driver_data, aes(x = time, 
                                          y = .data[[paste0(driver_name, "_detrendedNorm")]]), 
                  method = "loess", 
                  span = 0.1, 
                  se = FALSE, 
                  size = 1,
                  color = "black") +
      geom_line(data = abundance_data, aes(x = date, y = Anomaly_yr, 
                                           color = "Abundance"), size = 1.2) +
      annotate("text", x = as.POSIXct("1975-01-01"), y = 3.2,
               label = paste("r =", format(round(current_correlations[2], 2), 
                                           nsmall = 2)),
               parse = F, hjust = 0, vjust = 1, color = "black", size = 5.5) +
      annotate("text", x = as.POSIXct("1975-01-01"), y = 2.8,
               label = if (current_pvals[2] < 0.001) {
                 as.expression(bquote(italic(P)["conventional"] < 0.001))
               } else {
                 as.expression(bquote(italic(P)["conventional"] == .(round(current_pvals[2], 3))))
               },
               parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
      annotate("text", x = as.POSIXct("1975-01-01"), y = 2.4,
               label = if (current_boot_pvals[2] < 0.001) {
                 as.expression(bquote(italic(P)[bootstrap] < 0.001))
               } else {
                 as.expression(bquote(italic(P)[bootstrap] == .(round(current_boot_pvals[2], 3))))
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
    
    # Store the final plot for the current driver
    group_plots_list[[driver_name]] <- final_plot
  }
  # Store the final plots for the current group
  final_plots_list[[group_name]] <- group_plots_list
}
