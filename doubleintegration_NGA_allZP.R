################################################################################
#############          Pelagic Synthesis           #############################
#############             MAR-2025                 #############################
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
library(cowplot) #v1.1.3
library(here)
#library(tidyr)

## ------------------------------------------ ##
#            Biology Data -----
## ------------------------------------------ ##
ZP <- read.csv(file.path("raw",
                             "NGA",
                             "allzooplankton_NGA.csv")) %>%
  mutate(date = as.Date(paste0(Year, "-03-01")))  %>%
  arrange(taxa, season, Year) %>%
  group_by(taxa, season) %>%
  mutate(Anomaly_yr = approx(Year, Anomaly_yr, xout = Year)$y) %>%
  ungroup() %>%
  as.data.frame()

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

## ------------------------------------------ ##
#     Correlations -----
## ------------------------------------------ ##
# initialize lists
correlation_list <- list()
p_value_list <- list()

combinations <- unique(ZP[, c("season", "taxa")])

# loop over combinations of season and taxa
for (i in seq_len(nrow(combinations))) {
  season <- combinations[i, "season"]
  taxa <- combinations[i, "taxa"]
  
  cat("Processing group:", paste(season, taxa, sep = "_"), "\n")
  
  # filter data for current combo
  group_data <- ZP %>%
    filter(season == !!season, taxa == !!taxa)
  
  pdo_driver <- approx(pdoInt$time, pdoInt$pdoNorm, xout = group_data$date)$y
  interpolated_driver <- approx(pdoInt$time, pdoInt$pdoInt, xout = group_data$date)$y
  
  # calculate correlations
  correlations <- c(
    cor(group_data$Anomaly_yr, pdo_driver, use = "complete.obs"),
    cor(group_data$Anomaly_yr, interpolated_driver, use = "complete.obs")
  )
  # calculate pvalues 
  p_values <- c(
    cor.test(group_data$Anomaly_yr, pdo_driver)$p.value,
    cor.test(group_data$Anomaly_yr, interpolated_driver)$p.value
  )
  
  # save results
  label <- paste(season, taxa, sep = "_")
  correlation_list[[label]] <- correlations
  p_value_list[[label]] <- p_values
}

correlation_list
p_value_list

labels <- names(correlation_list)

# create df
results_df <- data.frame(
  label      = labels,
  cor_Norm   = sapply(correlation_list, `[`, 1),
  cor_Int    = sapply(correlation_list, `[`, 2),
  pval_Norm  = sapply(p_value_list,    `[`, 1),
  pval_Int   = sapply(p_value_list,    `[`, 2)
)

results_df <- results_df %>%
  separate(label, into = c("season", "taxa"), sep = "_")

rownames(results_df) <- NULL #remove row names
#write.csv(results_df, "output/NGA/CorAllZP_NGA.csv")

## ------------------------------------------ ##
#     PLOTS -----
## ------------------------------------------ ##
# list to store the final plots 
final_plots_list <- list()

output_dir <- here("figures/NGA/doubleint_allTaxa")

for (group_name in names(correlation_list)) {
  cat("Processing group:", group_name, "\n")
  
  parts <- strsplit(group_name, "_")[[1]]
  season <- parts[1]
  taxa <- parts[2]
  
  abundance_data <- ZP %>%
    filter(season == !!season, taxa == !!taxa)
  
  # pull correlation results
  cor_vals <- correlation_list[[group_name]]
  p_vals <- p_value_list[[group_name]]
  
  # legend label
  driver_legend <- "PDO"
  
  # --- Normalized Plot ---
  norm_plot <- ggplot() +
    # Driver time series (PDO normalized)
    geom_line(data = pdoInt, aes(x = time, y = pdoNorm, 
                                 color = "PDO"), linewidth = 1.2) +
    # Abundance time series
    geom_line(data = abundance_data, aes(x = date, y = Anomaly_yr, 
                                         color = "Abundance"), linewidth = 1.2) +
    # annotate correlation coefficient
    annotate("text", x = as.Date("2005-01-01"), y = 3,
             label = sprintf("r == %.2f", round(cor_vals[1], 2)),
             parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
    # annotate p-value
    annotate("text", x = as.Date("2005-01-01"), y = 2.4,
             label = if (p_vals[1] < 0.001) {
               as.expression(bquote(italic(P) < 0.001))
             } else {
               as.expression(bquote(italic(P) == .(round(p_vals[1], 3))))
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
    scale_x_date(breaks = seq(as.Date("1996-01-01"), 
                                  as.Date("2022-01-01"), 
                                  by = "4 years"),
                     date_labels = "%Y", expand = c(0, 0)) +
    coord_cartesian(xlim = as.Date(c("1996-01-01", "2022-01-01"))) +
    # manual color legend
    scale_color_manual(values = c("PDO" = "red", "Abundance" = "blue")) +
    labs(y = "PDO")
  
  # --- Integrated Plot ---
  int_plot <- ggplot() +
    # Driver time series (PDO integrated)
    geom_line(data = pdoInt, aes(x = time, y = pdoInt, 
                                 color = "PDO"), size = 1.2) +
    # Abundance time series
    geom_line(data = abundance_data, aes(x = date, y = Anomaly_yr, 
                                         color = "Abundance"), size = 1.2) +
    # annotate correlation coefficient
    annotate("text", x = as.Date("2005-01-01"), y = 3,
             label = sprintf("r == %.2f", round(cor_vals[2], 2)),
             parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
    # annotate p-value
    annotate("text", x = as.Date("2005-01-01"), y = 2.4,
             label = if (p_vals[2] < 0.001) {
               as.expression(bquote(italic(P) < 0.001))
             } else {
               as.expression(bquote(italic(P) == .(round(p_vals[2], 3))))
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
    scale_x_date(breaks = seq(as.Date("1996-01-01"), as.Date("2022-01-01"), 
                              by = "4 years"),
                 date_labels = "%Y", expand = c(0, 0)) +
    coord_cartesian(xlim = as.Date(c("1996-01-01", "2022-01-01"))) +
    # manual color legend
    labs(y = "Integrated PDO", color = "") +
    scale_color_manual(values = c("PDO" = "red", "Abundance" = "blue"))
  
  # extract and remove legend
  legend <- cowplot::get_legend(norm_plot)
  norm_plot <- norm_plot + theme(legend.position = "none")
  int_plot <- int_plot + theme(legend.position = "none")
  
  # title
  title <- ggdraw() +
    draw_label(label = bquote(paste(italic(.(taxa)), " in ", .(season))), 
               size = 18)
  
  # combine everything
  final_plot <- plot_grid(
    title,
    plot_grid(norm_plot, int_plot, ncol = 1, align = "v"),
    legend,
    ncol = 1, rel_heights = c(0.07, 0.85, 0.08)
  )
  
  print(final_plot)
  
  # save
  ggsave(
    filename = file.path(output_dir, paste0("PDO_NGA_", group_name, ".png")),
    plot = final_plot, width = 8, height = 8, dpi = 300, bg = "white"
  )
}
