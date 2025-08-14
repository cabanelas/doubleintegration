################################################################################
#############          Pelagic Synthesis           #############################
#############                 2025                 #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas 
################################################################################

# plots for ASLO ASM 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(tidyverse)
library(cowplot) #v1.1.3

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
final_results <- read.csv(file.path("output",
                                    "CCE",
                                    "boot_results_CCE.csv"), 
                                  header = T)

pval1 <- final_results$pval_pair[1]
pval2 <- final_results$pval_pair_integrated[1]

## correlations 
cor_CCE <- read.csv(file.path("output",
                              "CCE",
                              "weightedCor_CCE.csv"))
cor1 <- cor_CCE$cor_Norm[1]
cor2 <- cor_CCE$cor_Int[1]

## PDO 
PDO <- read.csv(file.path("output",
                          "CCE",
                          "PDO_CCE_integrated.csv"), 
                header = T) %>%
  select(-X) %>% 
  mutate(time = as.Date(time, 
                        format = "%Y-%m-%d"))

## Abundance
euphs <- read.csv(file.path("raw",
                            "CCE",
                            "nsimplex_CCE.csv")) %>%
  mutate(taxa = "Nsimplex") 

euphs$Yc_mean <- mean(euphs$Abundance)
euphs$Yc_sd <- sqrt(sum((euphs$Abundance - mean(euphs$Abundance, na.rm = TRUE))^2, na.rm = TRUE) / 
                      sum(!is.na(euphs$Abundance))) #pop mean; not sample..
#zscored
euphs$Anomaly_yr <- (euphs$Abundance - euphs$Yc_mean)/euphs$Yc_sd

#need to add the missing years 
#no sampling 1967, 1968, 1971, 1973 and 2020
#real 0s = 1972, 1976, 2010-2012
full_years <- data.frame(Year = seq(min(euphs$Year), max(euphs$Year), by = 1))
euphs <- full_years %>%
  left_join(euphs, by = "Year")
rm(full_years)

#linear interpolation 
#these analyses dont like NAs
euphs$Anomaly_yr <- approx(euphs$Year, 
                           euphs$Anomaly_yr,
                           xout = euphs$Year)$y

euphs <- euphs %>%
  fill(Region, Yc_mean, Yc_sd, taxa, .direction = "downup") 

# give abundance data day & month just to match w index
euphs$date <- as.Date(paste0(euphs$Year, "-03-01"))

## ------------------------------------------ ##
#     Plots -----
## ------------------------------------------ ##

## ------------------ ##
#     1) -----
## ------------------ ##
(NormPlot_bas <- ggplot() +
   # Driver time series (PDO normalized)
   geom_line(data = PDO, aes(x = time, y = pdoNorm, color = "PDO"), 
             linewidth = 1.2) +
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
         plot.title = element_blank(),
         panel.grid = element_blank()) +
   # x-axis settings
   scale_x_date(breaks = seq(as.Date("1950-01-01"), as.Date("2022-01-01"), 
                             by = "10 years"),
                date_labels = "%Y", expand = c(0, 0)) +
   coord_cartesian(xlim = as.Date(c("1950-01-01", "2022-01-01"))) +
   # manual color legend
   scale_color_manual(values = c("PDO" = "red", "Abundance" = "blue")) +
   labs(y = "PDO")
)


(Int1Plot_bas <- ggplot() +
    # Driver time series (Integrated PDO)
    geom_line(data = PDO, aes(x = time, y = pdoInt, color = "PDO"), 
              linewidth = 1.2) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 19, color = "black", face = "bold"),
          axis.text = element_text(size = 15, color = "black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 15, face = "bold"),
          legend.key.size = unit(2, "lines"),
          legend.position = "bottom",
          plot.title = element_blank(),
          panel.grid = element_blank()) +
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

legend_bas <- cowplot::get_plot_component(NormPlot_bas, 'guide-box-bottom', 
                                      return_all = TRUE)
cowplot::ggdraw(legend_bas)

# remove legend
NormPlot_bas <- NormPlot_bas + theme(legend.position = "none")
Int1Plot_bas <- Int1Plot_bas + theme(legend.position = "none")

# combine plots into one grid
combined_plots_bas <- plot_grid(NormPlot_bas, Int1Plot_bas,
                            ncol = 1, align = "v")

# add plot title 
title <- ggdraw() +
  draw_label(label = bquote("CCE " * italic("Nyctiphanes simplex") * " Spring"), 
             size = 18) +
  theme(plot.background = element_rect(fill = "transparent"))


# combine title and plots
(plot1_ASLO <- plot_grid(title, combined_plots_bas, legend_bas, ncol = 1, 
                         rel_heights = c(0.07, 0.8, 0.07)))

# export and save plot
ggsave("figures/CCE/bootstrap_CCE_Ns_spring_ASLO1.png", 
       plot1_ASLO, 
       width = 8, height = 8, dpi = 300, 
       bg = "white")



## ------------------ ##
#     2) -----
## ------------------ ##

(NormPlot <- ggplot() +
    # Driver time series (PDO normalized)
    geom_line(data = PDO, aes(x = time, y = pdoNorm, color = "PDO"), 
              linewidth = 1.2) +
    # Nsimplex time series (Abundance)
    geom_line(data = euphs, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
              linewidth = 1.2) +
    # annotate correlation coefficient
    annotate("text", x = as.Date("1960-01-01"), y = 3,
             label = sprintf("r == %.2f", cor1),
             parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
    # annotate bootstrapped p-value
    annotate("text", x = as.Date("1960-01-01"), y = 2.4,
             label = if (pval1 < 0.001) {
               as.expression(bquote(italic(p) < 0.001))
             } else {
               as.expression(bquote(italic(p) == .(round(pval1, 3))))
             },
             parse = TRUE, hjust = 0, vjust = 1, size = 5.5) +
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
          plot.title = element_blank(),
          panel.grid = element_blank()) +
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
    geom_line(data = PDO, aes(x = time, y = pdoInt, color = "PDO"), 
              linewidth = 1.2) +
    # Nsimplex time series (Abundance)
    geom_line(data = euphs, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
              linewidth = 1.2) +
    # annotate correlation coefficient
    annotate("text", x = as.Date("1960-01-01"), y = 3,
             label = sprintf("r == %.2f", cor2),
             parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
    # annotate bootstrapped p-value
    annotate("text", x = as.Date("1960-01-01"), y = 2.4,
             label = if (pval2 < 0.001) {
               as.expression(bquote(italic(p) < 0.001))
             } else {
               as.expression(bquote(italic(p) == .(round(pval2, 3))))
             },
             parse = TRUE, hjust = 0, vjust = 1, size = 5.5) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 19, color = "black", face = "bold"),
          axis.text = element_text(size = 15, color = "black"),
          legend.title = element_blank(),
          legend.text = element_text(size = 15, face = "bold"),
          legend.key.size = unit(2, "lines"),
          plot.title = element_blank(),
          panel.grid = element_blank()) +
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
ggsave("figures/CCE/bootstrap_CCE_Ns_spring.png", 
       final_plot, 
       width = 8, height = 8, dpi = 300, 
       bg = "white")








(NormPlot <- ggplot() +
    # Driver time series (PDO normalized)
    geom_line(data = PDO, aes(x = time, y = pdoNorm, color = "PDO"), 
              linewidth = 1.2) +
    # Nsimplex time series (Abundance)
    geom_line(data = euphs, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
              linewidth = 1.2) +
    # annotate correlation coefficient
    annotate("text", x = as.Date("1960-01-01"), y = 3,
             label = sprintf("rho == %.4f", cor1),
             parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
    # annotate bootstrapped p-value
    annotate("text", x = as.Date("1960-01-01"), y = 2.4,
             label = if (pval1 < 0.001) {
               as.expression(bquote(italic(P) < 0.001))
             } else {
               as.expression(bquote(italic(P) == .(round(pval1, 3))))
             },
             parse = TRUE, hjust = 0, vjust = 1, size = 5.5) +
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
    geom_line(data = PDO, aes(x = time, y = pdoInt, color = "PDO"), 
              linewidth = 1.2) +
    # Nsimplex time series (Abundance)
    geom_line(data = euphs, aes(x = date, y = Anomaly_yr, color = "Abundance"), 
              linewidth = 1.2) +
    # annotate correlation coefficient
    annotate("text", x = as.Date("1960-01-01"), y = 3,
             label = sprintf("rho == %.4f", cor2),
             parse = TRUE, hjust = 0, vjust = 1, color = "black", size = 5.5) +
    # annotate bootstrapped p-value
    annotate("text", x = as.Date("1960-01-01"), y = 2.4,
             label = if (pval2 < 0.001) {
               as.expression(bquote(italic(P) < 0.001))
             } else {
               as.expression(bquote(italic(P) == .(round(pval2, 3))))
             },
             parse = TRUE, hjust = 0, vjust = 1, size = 5.5) +
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
ggsave("figures/CCE/bootstrap_CCE_Ns_spring.png", 
       final_plot, 
       width = 8, height = 8, dpi = 300, 
       bg = "white")