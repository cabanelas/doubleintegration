################################################################################
#############          Pelagic Synthesis           #############################
#############                 2025                 #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas 
################################################################################
## NGA 
## Data exploration 

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(tidyverse)

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
files <- list.files("raw/NGA", 
                    pattern = "(_Fall|_Spr)\\.csv$", 
                    full.names = TRUE)

zp_all <- map_dfr(files, ~ {
  fname <- basename(.x)
  taxa <- str_extract(fname, "^[^B]+")
  season <- ifelse(str_detect(fname, "_Fall"), "Fall", "Spring")
  
  read_csv(.x, show_col_types = FALSE) %>%
    mutate(
      `Log STD` = suppressWarnings(as.numeric(`Log STD`)),
      CL95 = suppressWarnings(as.numeric(CL95)), #issues with some being read as character
      taxa = taxa, #add taxa column
      season = season, #add season column 
      #source_file = fname 
    )
})

## ------------------------------------------ ##
#            Tidy -----
## ------------------------------------------ ##
zp_all %>%
  filter(!is.na(`#samples`) & !is.na(`#collected`)) %>%
  summarise(all_equal = all(as.numeric(`#samples`) == as.numeric(`#collected`)))

# n and n-present same cols
# #samples and ...6 same cols
# #samples and #collected same cols
zp_all <- zp_all %>%
  mutate(
    n = coalesce(as.numeric(n), as.numeric(`n-present`)),
    `#samples` = coalesce(as.numeric(`#samples`), as.numeric(`...6`), as.numeric(`#collected`))
  ) %>%
  select(-c(`n-present`, `...6`, `#collected`)) %>%
  rename(num_samples = `#samples`,
         LogSTD = `Log STD`) %>%
# remove character from year column 
  mutate(Year = as.integer(str_extract(Year, "\\d{4}")))

# calculate anomalies zscore
zp_all <- zp_all %>%
  group_by(taxa, season) %>%
  mutate(
    Yc_mean = mean(LogMean, na.rm = TRUE),
    Yc_sd = sd(LogMean, na.rm = TRUE),
    Anomaly_yr = (LogMean - Yc_mean) / Yc_sd
  ) %>%
  ungroup()

#write.csv(zp_all, "raw/NGA/allzooplankton_NGA.csv")
## ------------------------------------------ ##
#            Plot TS -----
## ------------------------------------------ ##
for (taxon in unique(zp_all$taxa)) {

  p <- zp_all %>%
    filter(taxa == taxon) %>%
    ggplot(aes(x = Year, y = Anomaly_yr)) + 
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_line() +
    geom_point() +
    facet_wrap(~season, ncol = 1) +
    labs(
      title = paste(taxon),
      y = "Anomaly (Z-score)",
      x = "Year"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 15, color = "black"))
  
  print(p)
  
  #ggsave(filename = paste0("figures/NGA/", 
  #                         taxon, "_anomaly_timeseries.png"),
  #       plot = p, width = 6, height = 5, 
  #       bg = "white")
}
