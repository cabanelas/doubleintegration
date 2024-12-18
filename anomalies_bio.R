################################################################################
#############          Pelagic Synthesis           #############################
#############             MAR-2024                 #############################
#############     Calculating ECOMON Anomalies     #############################
## by: Alexandra Cabanelas 
################################################################################
## Double Integration Analysis NES
# Script #1 
# script to calculate bio anomalies
#     a. Log(x + min/2) transformation
#     b. Spatial average
#     c. Z-score transformation

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##

library(tidyverse) #v2.0.0
library(listviewer) #v4.0.0; for looking at lists interactively 
#library(magrittr) #v2.0.3; map_dfr
library(here) #v1.0.1; easily build path to files  

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##

# can download EcoMon data from 
#https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0187513 
# NCEI Accession 0187513 v3.3 

abu <- read.csv(file.path("raw","EcoMon_v3_8_wDateStrata.csv"))

# preconditions: raw EcoMon data was processed to georeference based on the two 
# provided shapefiles: EcomonStrata_v4.shp & EcomonStrata_v4b.shp
# to classify the distinct ecoregions within the NES using the provided 108 
# spatial polygons (e.g., identify Gulf of Maine, Georges Bank, 
# Southern New England, and Mid-Atlantic Bight)

## ------------------------------------------ ##
#            Tidy Data -----
## ------------------------------------------ ##

# using 10m2 values; [100m3 data is available]
names(abu)<-gsub("_10m2","",names(abu)) # get rid of "_10m2" in colnames

# Select taxa to analyze
#c.typicus, c.finmarchicus, pseudocalanus
taxa_of_interest <- c("ctyp", "calfin", "pseudo")

# add season to df
abu <- abu %>%
  mutate(season = case_when(between(month, 3, 5) ~ "spring",
                            between(month, 6, 8) ~ "summer",
                            between(month, 9, 11) ~ "fall",
                            TRUE ~ "winter"))
# add region to df
abu <- abu %>%
  mutate(Region = case_when(region == 1 ~ "MAB", #MidAtlantic Bight
                            region == 2 ~ "SNE", #Southern New England
                            region == 3 ~ "GB", #Georges Bank
                            region == 4 ~ "GOM", #Gulf of Maine
                            TRUE ~ "Outside"))

# select cols of interest
abu1 <- abu %>%
  select(date, month, day, year, all_of(taxa_of_interest), season, Region)

# pivot from wide to long
abu_long <- abu1 %>%
  pivot_longer(cols = c(ctyp:pseudo), #adjust as needed, depending on taxa used
               names_to = "taxa", values_to = "abundance")

#remove rows with NANS = no sampling for zp or itchyo 
#zoo_not_ich = 3437
#ich_not_zoo = 3627
#both = 25629          
#both gears not always used
abu_long <- abu_long %>% filter(!is.nan(abundance))

################################################################################
# list of taxa and regions - for the loop below
# not essential; used it for regions to exclude "outside" 
# can use it to select specific taxa; but we also filtered those out above in 
#line 47

#taxa_list <- c("ctyp", "calfin")  
region_list <- c("MAB", "SNE", "GOM", "GB")  #excluding "outside"
#season_list <- c("spring","summer","winter","fall")

################################################################################
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
################################################################################
################################################################################
# 1) Biological TS
#     a. Log(x + min/2) transformation
#     b. Spatial average
#     c. Z-score transformation for plotting on same axes = anomalies 

## ------------------------------------------ ##
#     Calculate anomalies/biological TS -----
## ------------------------------------------ ##

zscore_list <- list()

for (taxa in unique(abu_long$taxa)) {
  for (Region in region_list) { #excluding "outside" region 
    for (season in unique(abu_long$season)) {
      cat("Processing:", taxa, "-", Region, "-", season, "\n")
      
      abu_long_subset <- abu_long[abu_long$taxa == taxa & 
                                    abu_long$Region == Region & 
                                    abu_long$season == season, ]
      
      #print(dim(abu_long_subset))
      
      # Step 1a: Data transform
      abu_longTran <- abu_long_subset %>%
        group_by(taxa, Region, season) %>%
        # Find the min. non-zero value for data transform for each region and season
        mutate(min_nonzero = min(abundance[abundance != 0], na.rm = TRUE)) %>%
        ungroup() %>%
        #to find the log for each year, have to group by yr
        group_by(taxa, Region, season, year) %>%
        # Log transformation - common/base10
        mutate(LogM = log10(abundance + min_nonzero/2)) %>%
        ungroup()
      
      # Step 1b: Spatial average
      # 1b.1: Calculate the mean of LogM values by season, year, taxa, and region
      abu_longMean <- abu_longTran %>%
        group_by(season, year, taxa, Region) %>%
        summarize(mean_LogM = mean(LogM, na.rm = TRUE))
      
      # 1b.2: Calculate the standard deviation of LogM values by season, taxa, and region
      abu_longSD  <- abu_longMean %>%
        group_by(season, taxa, Region) %>% #not year, to get whole ts 
        summarize(Yc_sd1 = sd(mean_LogM, na.rm = TRUE),
                  Yc_mean1 = mean(mean_LogM, na.rm = TRUE)) 
      
      # Merge with abu_longMean to include year information in final df
      abu_longSD <- abu_longSD %>%
        left_join(abu_longMean, by = c("season", "taxa", "Region")) 
      
      # Step 1c: Calculate z-score/anomalies
      zscore1 <- abu_longSD %>%
        mutate(Anomaly_yr = (mean_LogM - Yc_mean1) / Yc_sd1)
      
      # Merge with abu_longTran to include the date column
      zscore2 <- zscore1 %>%
        left_join(abu_longTran %>% select(season, year, date), 
                  by = c("season", "year"))
      
      # Convert date column to POSIXct
      zscore2$date <- as.Date(zscore2$date, format = "%m/%d/%Y")
      
      zscore <- zscore2 %>%
        distinct(season,taxa,Region,year, .keep_all = T)
      
      # Store zscore in the list
      zscore_list[[paste(taxa, Region, season, sep = "_")]] <- zscore
    }
  }
}

zscore <- do.call(rbind, zscore_list) 
#write.csv(zscore, "output/zscore_anomalies_3taxazp_28MAY.csv")
rm(zscore1, zscore2, abu_longMean, abu_longSD, abu_longTran)



ggplot(zscore, aes(x=year, y=Anomaly_yr)) + 
  geom_line() + 
  facet_grid(taxa ~ season+Region)

#zscore_list <- zscore_list[-which(names(zscore_list) == "ammspp_SNE_fall")]

# check SD and mean of anomalies by group
zscore %>%
  group_by(taxa, Region, season) %>%
  summarize(mean_Anomaly_yr = mean(Anomaly_yr, na.rm = TRUE),
            sd_Anomaly_yr = sd(Anomaly_yr, na.rm = TRUE)) %>%
  print(n = 50)

# remove ammspp fall SNE since == 0 
zscore <- zscore[complete.cases(zscore$Anomaly_yr), ]

#checking the list output from the loop - making sure it looks right
listviewer::jsonedit(zscore_list)
purrr::map(zscore_list, "name") #to see taxa_region_season
purrr::map(zscore_list, 2) #to check that taxa were properly assigned 
purrr::map(zscore_list, 1) #can also pipe this 

# to extract as tibble
#map_dfr(zscore_list, extract, c("season", "taxa", "Region", "Yc_sd1",
                                #"Yc_mean1","year","mean_LogM","Anomaly_yr","date"))
#zscore_list %>% pluck("ctyp_MAB_spring") #to access specific 'groups'
