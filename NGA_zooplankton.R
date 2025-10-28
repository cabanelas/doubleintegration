################################################################################
#############          Pelagic Synthesis           #############################
#############   NGA - Double Integration Analysis  #############################
## by: Alexandra Cabanelas 
## created MAR-2025, updated OCT-2025

################################################################################
## Northern Gulf of Alaska LTER
## joining all separate taxa csv into one
## Fall 1997-2022 & Spring 1998-2022
# small taxa = m3; large taxa = m2
# "n-present" = # of stations that had data and are used to calculate mean and stats
# "#collected = # of stations samples (includes the 0s)

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(tidyverse)
library(zoo)

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
files <- list.files("raw/NGA", 
                    pattern = "(_Fall|_Spr|_Sprv2)\\.csv$", 
                    full.names = TRUE)

zp_all <- map_dfr(files, ~ {
  fname <- basename(.x)
  # treat both _Spr and _Sprv2 as Spring
  # have to do this bc I received updated file for NeocalanusCris_Spr
  season <- ifelse(str_detect(fname, "_Fall"), "Fall", "Spring")
  # extract taxa name from file name
  taxa <- str_extract(fname, "^[^B]+")
  read_csv(.x, show_col_types = FALSE) %>%
    mutate(
      `Log STD` = suppressWarnings(as.numeric(`Log STD`)),
      CL95 = suppressWarnings(as.numeric(CL95)), #issues with some being read as character
      taxa = taxa,
      season = season
    )
})

## ------------------------------------------ ##
#            Tidy -----
## ------------------------------------------ ##
# check whether values in the "#samples" col are == values in "#collected"
zp_all %>%
  filter(!is.na(`#samples`) & !is.na(`#collected`)) %>%
  summarise(all_equal = all(as.numeric(`#samples`) == as.numeric(`#collected`)))

# n and n-present same cols
# #samples and ...6 same cols
# #samples and #collected same cols
zp_all <- zp_all %>%
  mutate(
    # if n is missing, fill it with the value from n-present
    n = coalesce(as.numeric(n), as.numeric(`n-present`)),
    # if #samples is missing, fill it with the value from ...6 or #collected
    `#samples` = coalesce(as.numeric(`#samples`), as.numeric(`...6`), as.numeric(`#collected`))
  ) %>%
  select(-c(`n-present`, `...6`, `#collected`)) %>%
  rename(num_samples = `#samples`,
         LogSTD = `Log STD`) %>%
  # remove character from year column 
  mutate(Year = as.integer(str_extract(Year, "\\d{4}")))

## --- check if there are any missing years to interpolate ------ 
# define year range for df by season
# Fall 1997-2022 & Spring 1998-2022
# larvacean spring + fall 2001-2021
# limacina spring + fall 2001-2021
check_missing <- zp_all %>%
  group_by(taxa, season) %>%
  summarise(
    years_present = list(sort(unique(Year))),
    n_na = sum(is.na(LogMean)),
    .groups = "drop"
  ) %>%
  mutate(
    expected_years = case_when(
      taxa %in% c("Larvacean", "Limacina") ~ list(2001:2021),
      season == "Fall"   ~ list(1997:2022),
      season == "Spring" ~ list(1998:2022),
      TRUE ~ list(integer(0))
    ),
    missing_years = map2(expected_years, years_present, setdiff),
    n_missing = map_int(missing_years, length)
  ) %>%
  filter(n_missing > 0 | n_na > 0)

## --- linear interpolate the missing LogMean ------
zp_all1 <- zp_all %>%
  group_by(taxa, season) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(LogMean = zoo::na.approx(LogMean, Year, na.rm = FALSE)) %>%
  ungroup()

# check which group was interpolated
zp_all1 %>%
  left_join(zp_all, by = c("taxa", "season", "Year"), 
            suffix = c("_filled", "_orig")) %>%
  filter(is.na(LogMean_orig) & !is.na(LogMean_filled)) %>%
  select(taxa, season, Year, LogMean_filled)

## ------------------------------------------ ##
#       Calculate anomalies zscore -----
## ------------------------------------------ ##
zp_all1 <- zp_all1 %>%
  group_by(taxa, season) %>%
  mutate(
    Yc_mean = mean(LogMean, na.rm = TRUE),
    Yc_sd = sd(LogMean, na.rm = TRUE),
    Anomaly_yr = (LogMean - Yc_mean) / Yc_sd
  ) %>%
  ungroup()

#write.csv(zp_all1, "raw/NGA/allzooplankton_NGA.csv")
## ------------------------------------------ ##
#            Plot TS -----
## ------------------------------------------ ##
for (taxon in unique(zp_all1$taxa)) {

  p <- zp_all1 %>%
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
