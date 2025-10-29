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
library(here)
#library(tidyr)

## ------------------------------------------ ##
#           1) Biology Data -----
## ------------------------------------------ ##
ZP <- read_csv(file.path("raw",
                         "NGA",
                         "allzooplankton_NGA.csv")) %>%
  select(, -1) %>%
  mutate(date = as.Date(paste0(Year, "-03-01"))) #%>% #give abundance daymonth to match w index
  #arrange(taxa, season, Year) %>%
  #group_by(taxa, season) %>%
  #mutate(Anomaly_yr = approx(Year, Anomaly_yr, xout = Year)$y) %>%
  #ungroup() 

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

# initialize lists
correlation_list <- list()
p_value_list <- list()

combinations <- unique(ZP[, c("season", "taxa")])

# loop over combinations of season and taxa
for (i in seq_len(nrow(combinations))) {
  #season <- combinations[i, "season"]
  #taxa <- combinations[i, "taxa"]
  season <- combinations$season[i]  
  taxa   <- combinations$taxa[i]    
  
  cat("Processing group:", paste(season, taxa, sep = "_"), "\n")
  
  # filter data for current combo
  group_data <- ZP %>%
    dplyr::filter(season == !!season, taxa == !!taxa)
  
  pdo_driver <- approx(PDO$time, PDO$pdo_z, xout = group_data$date)$y
  interpolated_driver <- approx(PDO$time, PDO$pdo_int_z, xout = group_data$date)$y
  
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
  site = "NGA",
  label      = labels,
  cor_Norm   = sapply(correlation_list, `[`, 1),
  cor_Int    = sapply(correlation_list, `[`, 2),
  pval_Norm  = sapply(p_value_list,    `[`, 1),
  pval_Int   = sapply(p_value_list,    `[`, 2),
  tau_months = tau_bio
)

results_df <- results_df %>%
  separate(label, into = c("season", "taxa"), sep = "_")

rownames(results_df) <- NULL #remove row names
#write.csv(results_df, "output/NGA/cor_All_ZP_NGA.csv")

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
    geom_line(data = PDO, aes(x = time, y = pdo_z, 
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
               as.expression(bquote(italic(p) < 0.001))
             } else {
               as.expression(bquote(italic(p) == .(round(p_vals[1], 3))))
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
    scale_x_date(breaks = seq(as.Date("1995-01-01"), 
                                  as.Date("2022-01-01"), 
                                  by = "5 years"),
                     date_labels = "%Y", expand = c(0, 0)) +
    coord_cartesian(xlim = as.Date(c("1995-01-01", "2022-01-01"))) +
    # manual color legend
    scale_color_manual(values = c("PDO" = "red", "Abundance" = "blue")) +
    labs(y = "PDO")
  
  # --- Integrated Plot ---
  int_plot <- ggplot() +
    # Driver time series (PDO integrated)
    geom_line(data = PDO, aes(x = time, y = pdo_int_z, 
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
               as.expression(bquote(italic(p) < 0.001))
             } else {
               as.expression(bquote(italic(p) == .(round(p_vals[2], 3))))
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
    scale_x_date(breaks = seq(as.Date("1995-01-01"), as.Date("2022-01-01"), 
                              by = "5 years"),
                 date_labels = "%Y", expand = c(0, 0)) +
    coord_cartesian(xlim = as.Date(c("1995-01-01", "2022-01-01"))) +
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
  #ggsave(
  #  filename = file.path(output_dir, paste0("PDO_NGA_", group_name, ".png")),
  #  plot = final_plot, width = 8, height = 8, dpi = 300, bg = "white"
  #)
}
