################################################################################
#############        LTER Pelagic Synthesis WG     #############################
#############        Double Integration Analysis   #############################
#############        Driver AR(1) coefficient      #############################
## by: Alexandra Cabanelas 
## created OCT-2025
################################################################################
### --- Drivers
#https://psl.noaa.gov/enso/dashboard.html
#https://psl.noaa.gov/data/timeseries/month/

## California Current Ecosystem LTER
# Pacific Decadal Oscillation          = PDO   == 1854 - 2024
# North Pacific Gyre Oscillation       = NPGO  == 1950 - 2025 https://o3d.org/npgo/

## Northern Gulf of Alaska LTER
# Oceanic Nino Index                   = ONI   == 1950 - DEC 2024
# Gulf of Alaska Downwelling Index     = GOADI == 1993 - 2024
# Northern Gulf of Alaska Oscillation  = NGAO  == 1993 - 2024
# Pacific North-American Pattern       = PNA   == 1950 - SEP 2023

## Palmer LTER
# Multivariate ENSO Index              = MEI   == 1983 - 2024
# Southern Annular Mode                = SAM   == 1957 - AUG 2023
## -------------------------------------------------------------------------- ##
# Script #1 : 01_ARcoef_ALL_drivers
# script to calculate AR coefficient of drivers - loop for all drivers 

## STEP 1 of Monte Carlo analysis
# Step 1 - estimate autoregression coeff from the original data/signals to be
#able to use for creating two red-noise time series 
# the goal is to use the AR1 coefficient to create red-noise surrogates that mimic
# the autocorreltion structure of the original PDO

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(forecast) #v8.21; Arima()
library(astsa) #v2.1; acf2 (optional); #library(urca)
library(tidyverse)
library(lubridate)

## ------------------------------------------ ##
#            Data & Tidy-----
## ------------------------------------------ ##

### --- California Current Ecosystem LTER -----
## --- PDO -----
PDO <- read_csv(file.path("raw", 
                          "CCE",
                          "PDO.csv")) %>% # contains up to DEC-2024
  # long format
  pivot_longer(cols = Jan:Dec, 
               names_to = "month", 
               values_to = "pdo") %>%
  mutate(
    # turn month into integer
    month = match(month, month.abb),
    # make date format
    date = make_date(year = Year, month = month, day = 1)
  ) %>%
  arrange(date) %>%
  # filtering 10 yrs before bio data (doesn't make much diff for ARcoef)
  filter(Year >= 1941, Year <= 2021, pdo < 99) %>%
  set_names(tolower)

## --- NPGO -----
# missing value -9999 https://psl.noaa.gov/data/timeseries/month/
NPGO <- read_csv(file.path("raw", 
                           "CCE",
                           "NPGO.csv")) %>% # contains up to MAR-2025
  separate(Date, 
           into = c("month", "day", "year"), 
           sep = "-", remove = FALSE) %>%
  mutate(
    # make date format
    year = as.integer(year),
    # force 2-digit years into correct century
    year = if_else(year >= 50, 1900 + year, 2000 + year),  
    Date = make_date(year, as.integer(month), as.integer(day))
  ) %>%
  # filtering 10 yrs before bio data (doesn't make much diff for ARcoef)
  filter(year <= 2021, npgo < 99) %>%
  set_names(tolower) %>%
  select(-day)

## --- ENSO ----- NEED TO ADD ENSO STUFF!!!! 
## ------------------------------------------ ##

### --- Northern Gulf of Alaska LTER -----
## --- ONI -----
ONI <- read_csv(file.path("raw",
                          "NGA",
                          "driver",
                          "ONI.csv")) %>% # contains up to DEC-2024
  rename(year = YR, oni = ANOM) %>% 
  select(-TOTAL) %>%
  mutate(month = case_when(
    SEAS == "DJF" ~ 1,
    SEAS == "JFM" ~ 2,
    SEAS == "FMA" ~ 3,
    SEAS == "MAM" ~ 4,
    SEAS == "AMJ" ~ 5,
    SEAS == "MJJ" ~ 6,
    SEAS == "JJA" ~ 7,
    SEAS == "JAS" ~ 8,
    SEAS == "ASO" ~ 9,
    SEAS == "SON" ~ 10,
    SEAS == "OND" ~ 11,
    SEAS == "NDJ" ~ 12
  ),
  # fix date
  date = make_date(year = year, month = month, day = 1) 
  ) %>%
  #filtering 10 yrs before bio data (doesn't make a diff for ARcoef)
  filter(year > 1987) %>%
  set_names(tolower)

## --- GOADI -----
GOADI <- read_csv(file.path("raw",
                            "NGA",
                            "driver",
                            "GOADI_monthly.csv")) %>% # contains up to SEP-2024
  select(-1) %>%
  rename(goadi = DW) %>% 
  # fix date
  mutate(
    Date = ym(Date), # convert to date object 
    Year = year(Date),          
    Month = month(Date)          
  ) %>%
  set_names(tolower)

## --- NGAO -----
NGAO <- read_csv(file.path("raw",
                           "NGA",
                           "driver",
                           "NGAO_monthly.csv")) %>% # contains up to SEP-2024
  select(-1) %>%
  # fix date
  mutate(
    Date = ym(Date), # convert to date object 
    Year = year(Date),          
    Month = month(Date)          
  ) %>%
  set_names(tolower)

## --- PNA -----
PNA <- read_csv(file.path("raw",
                          "NGA",
                          "driver",
                          "PNA.csv")) %>% # contains up to SEP-2023
  pivot_longer(cols = Jan:Dec, 
               names_to = "month", 
               values_to = "pna") %>%
  mutate(
    # turn month into integer
    month = match(month, month.abb),
    # make date format
    Date = make_date(year = year, month = month, day = 1)
  ) %>%
  drop_na(pna) %>% #drop rows where pna is NA
  arrange(Date) %>%
  set_names(tolower)
## ------------------------------------------ ##

### --- Palmer LTER -----
## --- MEI -----
MEI <- read.csv(file.path("raw", 
                          "PAL",
                          "MEI.csv")) %>% # contains up to FEB-2024
  mutate(
    DATE = mdy(DATE),   # make date format
    Year = year(DATE),  # extract year
    Month = month(DATE) # extract month
  ) %>%
  filter(Year > 1982) %>% #filtering 10 yrs before bio data
  set_names(tolower)

## --- SAM -----
SAM <- read_csv(file.path("raw", 
                          "PAL",
                          "SAM.csv")) %>% # contains up to aug-2023
  set_names(~ str_to_title(tolower(.x))) %>%
  pivot_longer(cols = Jan:Dec, 
               names_to = "month", 
               values_to = "sam") %>%
  mutate(
    # turn month into integer
    month = match(month, month.abb),
    # make date format
    Date = make_date(year = Year, month = month, day = 1)
   ) %>%
  filter(Year > 1982) %>% #filtering 10 yrs before bio data
  drop_na(sam) %>%
  set_names(tolower)
## ------------------------------------------ ##

## ------------------------------------------ ##
#            Loop -----
## ------------------------------------------ ##
PDO   <- PDO   %>% rename(value = pdo)
NPGO  <- NPGO  %>% rename(value = npgo)
ONI   <- ONI   %>% rename(value = oni)
GOADI <- GOADI %>% rename(value = goadi)
NGAO  <- NGAO  %>% rename(value = ngao)
PNA   <- PNA   %>% rename(value = pna)
MEI   <- MEI   %>% rename(value = mei)
SAM   <- SAM   %>% rename(value = sam)

drivers_list <- list(
  ONI   = ONI,
  GOADI = GOADI,
  NGAO  = NGAO,
  PDO   = PDO,
  NPGO  = NPGO,
  PNA   = PNA,
  MEI   = MEI,
  SAM   = SAM
)

# check colnames for all 
walk2(names(drivers_list), drivers_list, function(name, df) {
  cat("\n---", name, "---\n")
  print(colnames(df))
})

drivers_cleaned <- tribble(
  ~driver, ~df_object, ~start_year, ~start_month,
  "ONI",   "ONI",      1988,         1,
  "GOADI", "GOADI",    1993,         1,
  "NGAO",  "NGAO",     1993,         1,
  "PDO",   "PDO",      1941,         1,
  "NPGO",  "NPGO",     1950,         1,
  "PNA",   "PNA",      1988,         1,
  "MEI",   "MEI",      1983,         1,
  "SAM",   "SAM",      1983,         1
)

ar_info_df <- drivers_cleaned %>%
  pmap_df(function(driver, df_object, start_year, start_month) {
    df <- get(df_object)
    
    # --- create time series object -----
    ts_data <- ts(df$value, 
                  start = c(start_year, start_month), 
                  frequency = 12)
    
    # --- fit AR1 model -----
    AR1_model <- Arima(ts_data, order = c(1, 0, 0))
    
    # --- extract diagnostics -----
    tibble(
      driver     = driver,
      AR_coef    = AR1_model$coef["ar1"],
      se         = sqrt(diag(vcov(AR1_model)))[1],
      sigma2     = AR1_model$sigma2,
      n          = length(ts_data),
      year_start = start_year,
      year_end   = max(df$year)
    )
  })
#write.csv(ar_info_df, "output/ARcoef_ALL_drivers.csv")
## ------------------------------------------ ##


## ------------------------------------------ ##
#     diagnostic plots -----
## ------------------------------------------ ##

walk2(names(drivers_list), drivers_list, function(name, df) {
  ts_data <- ts(df$value, start = c(min(df$year), min(df$month)), frequency = 12)
  AR1_model <- Arima(ts_data, order = c(1, 0, 0))
  residuals_ar1 <- residuals(AR1_model)
  
  cat("\n---", name, "---\n")
  
  # plot diagnostics
  par(mfrow = c(2, 2))  
  
  # time series with fitted values
  ts.plot(ts_data, main = paste(name, "Time Series with AR(1) Fit"), ylab = name)
  lines(fitted(AR1_model), col = "red", lwd = 2)
  
  # residuals
  plot(residuals_ar1, main = "Residuals of AR(1) Model", ylab = "Residuals", type = "l")
  
  # histogram
  hist(residuals_ar1, main = "Histogram of Residuals", xlab = "Residuals")
  
  # Q-Q plot
  qqnorm(residuals_ar1)
  qqline(residuals_ar1, col = "red")
  
  # ACF and PACF
  acf2(residuals_ar1)
  par(mfrow = c(1, 1))
})

