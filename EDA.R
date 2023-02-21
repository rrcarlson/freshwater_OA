######### Data exploration
######### 07 February 2023
######### Rachel R. Carlson

# Libraries
library(tidyverse)
library(ncdf4)
library(raster)
library(sf)
library(purrr)
library(lubridate)

### Sorting datasets for examining freshwater impact on OA in California, Washington, Oregon.
# Filter by datasets with time series inside relevant dates
# full is loaded from RData `clean_all_distance_offshore_seacarb`
# Join datasets
full <- clean_all
OOI <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/66_Final_OOI_SLHNov2022.csv") %>% dplyr::select(-X)# OOI data (dataset 66)
Trinidad <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/67_TrinidadHeadLine/67_TrinidadHead_Nov2022SLH.csv") %>% dplyr::select(-X)
Newport <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/68_NewportHydrographicLine/68_NewportHydrographicLine_Final.csv") %>% dplyr::select(-c(X,X.1))
Columbia <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/70_CMOPSaturn2/70_CMOPSaturn2_final.csv") %>% dplyr::select(-X)
Scripps <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/71_ScrippsLaJolla/71_ScrippsLaJolla_final.csv") %>% dplyr::select(-X)
variable.names(full)
variable.names(Scripps)

# Dataset 27 pH values are filtered out of the 'clean' dataset because they use YSI. However, this dataset will be important to us.
# Filter out "clean" version of dataset 27 and then load the full (unfiltered) version back in.
full <- full %>% filter(dataset_id != 27)
CoastOA <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/27_CoastalOA/27_final_coastal_OA.csv")

# Bind datasets
OOI$distance_offshore <- NA # Create a blank column for distance_offshore as this is missing in individual datasets
Trinidad$distance_offshore <- NA
Newport$distance_offshore <- NA
Columbia$distance_offshore <- NA
Scripps$distance_offshore <- NA
CoastOA$distance_offshore <- NA

full2 <- do.call("rbind", list(full, OOI, Trinidad, Newport, Columbia, Scripps, CoastOA)) # Bind datasets

# id_list <- c(30, 31, 32, 33, 39, 40, 50, 56, 64, 66, 67, 68, 70, 71)
# trim <- full2 %>% filter(dataset_id %in% id_list)
# trim <- trim %>% group_by(dataset_id) %>% mutate(min_time = min(time_utc, na.rm = TRUE),
#                                             max_time = max(time_utc, na.rm = TRUE))
# 
# View(trim %>% group_by(dataset_id) %>% summarize(min_time = min(time_utc, na.rm = TRUE),
#                                          max_time = max(time_utc, na.rm = TRUE)))

# Visualize NA values.
# View(full2[is.na(full2$longitude),]) # It appears they are all in dataset 67, which is the Trinidad headlands. Monitoring sites are not near freshwater and could be complicated by upwelling, so we'll remove this dataset.
# View(full2[is.na(full2$latitude),])
# trim <- trim %>% filter(dataset_id != 67)
# # Map the sites
# coords1 <- trim %>% distinct(latitude, longitude, dataset_id, min_time, max_time)
# # Convert to shapefile
# mcoords <- st_as_sf(coords1, coords = c("longitude", "latitude"), crs = 4326)
# plot(mcoords$geometry)
# 
# st_write(mcoords,"/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/recommended_sites.9Feb2023.shp")

# The above (only those datasets that were recommended specifically for this study) seems to be missing a lot - Recalculate time ranges and and plot full data
full2 <- full2 %>% group_by(dataset_id) %>% mutate(min_time = min(time_utc, na.rm = TRUE),
                                                   max_time = max(time_utc, na.rm = TRUE))
coords_full <- full2 %>% distinct(latitude, longitude, dataset_id, min_time, max_time) %>% na.omit()
# Filter for datasets that are over 3 months (90 days) in time range
coords_full$time_diff <- coords_full$max_time - coords_full$min_time # Difference between youngest and oldest time
coords_ts <- coords_full %>% filter(time_diff > 90)
# Remove CalCOFI (dataset 25) because not nearshore
coords_ts <- coords_ts %>% filter(dataset_id != 25)
# Clean up errors in coordinates. There are some crazy values with inconsistent spatial extents that alter final mapping
coords_ts <- coords_ts %>% filter(latitude < 50)
coords_ts <- coords_ts %>% filter(longitude > -126)
  
mcoords_ts <- st_as_sf(coords_ts, coords = c("longitude", "latitude"), crs = 4326)

st_write(mcoords_ts,"/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/all_sites.9Feb2023.shp")

# I used a KML of the above to tag all sites located near river mouths. Need to filter sf object to this subset of sites.
# Will use lat.lon to filter data so first need to remove tailing zeros from lat.lon decimals
full2$latitude <- as.numeric(sub("0+$", "", as.character(full2$latitude))) %>% round(digits = 4)
full2$longitude <- as.numeric(sub("0+$", "", as.character(full2$longitude))) %>% round(digits = 4)

# Import sites tagged for river proximity
fresh <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/sites.csv")
fresh$Lat <- as.numeric(sub("0+$", "", as.character(fresh$Lat))) %>% round(digits = 4)
fresh$Lon <- as.numeric(sub("0+$", "", as.character(fresh$Lon))) %>% round(digits = 4)
fresh <- fresh %>% dplyr::select(c("Site", "Dataset_ID","Distance","Reach_ID","Flow_range_midFeb","Mean_flow","Lat","Lon"))
full2 <- full2 %>% rename("Lat" = "latitude", "Lon" = "longitude", "Dataset_ID" = "dataset_id") # There were 70ish mismatched dataset IDs in first try at joining these datasets. Ensure datest ID in tagged freshwater dataset is correct.

# Join reach, site name, flow, and distance information to relevant sites
full_fresh <- full2 %>% right_join(fresh, by = c("Lat" = "Lat", "Lon" = "Lon", "Dataset_ID" = "Dataset_ID"), multiple = "all")

# Filter by sites with pH measured
full_fresh_pH <- full_fresh %>% filter(!is.na(pH_total))
full_fresh_alk <- full_fresh %>% filter(!is.na(ta_umolkg))

# Check whether these sites have alkalinity, temp, salinity as well
full_fresh_4 <- full_fresh_pH %>% 
  filter(!is.na(sal_pss)) %>% 
  filter(!is.na(t_C)) %>% 
  filter(!is.na(ta_umolkg))

full_fresh_3 <- full_fresh_pH %>% 
  filter(!is.na(sal_pss)) %>% 
  filter(!is.na(t_C))

# Filter by sites with complete values for temp and pH to check for upwelling signal (or departures thereof)
full_fresh_pH_temp <- full_fresh_pH %>% 
  filter(!is.na(t_C))

# Tally number of observations per Site with temp and pH data
View(full_fresh_pH_temp %>% group_by(Site) %>% tally())
View(full_fresh_4 %>% group_by(Site) %>% tally())
View(full_fresh_3 %>% group_by(Site) %>% tally())

# Plot pH v. temp. Departures should indicate areas where upwelling is not explanatory
# Tomales Bay
# NO clear pattern (upwelling)
full_fresh %>% 
  filter(Site == "Tomales Bay") %>% 
  ggplot(aes(x = t_C, y = pH_total)) +
  geom_point()

# strong freshwater - alkalinity pattern
full_fresh %>% 
  filter(Site == "Tomales Bay") %>% 
  ggplot(aes(x = sal_pss, y = ta_umolkg)) +
  geom_point()

# Santa Barbara Pier - moderate flow
# NO clear pattern
full_fresh %>% 
  filter(Site == "Santa Barbara Pier; Mission Creek and Sycamore Creek") %>% 
  ggplot(aes(x = t_C, y = pH_total)) +
  geom_point()

# some freshwater - alkalinity pattern
full_fresh %>% 
  filter(Site == "Santa Barbara Pier; Mission Creek and Sycamore Creek") %>% 
  ggplot(aes(x = sal_pss, y = ta_umolkg)) +
  geom_point()

# Refugio - very low flow
# SMALL pattern from upwelling
full_fresh %>% 
  filter(Site == "Refugio") %>% 
  ggplot(aes(x = t_C, y = pH_total)) +
  geom_point()

# Moderate freshwater - alkalinity pattern
full_fresh %>% 
  filter(Site == "Refugio") %>% 
  ggplot(aes(x = sal_pss, y = ta_umolkg)) +
  geom_point()

# Morro Bay - medium flow
# NO distinct pattern
full_fresh %>% 
  filter(Site == "Morro Bay (Chorro Creek and Los Osos Creek)") %>% 
  ggplot(aes(x = t_C, y = pH_total)) +
  geom_point()

# alkalinity data not there
full_fresh %>% 
  filter(Site == "Morro Bay (Chorro Creek and Los Osos Creek)") %>% 
  ggplot(aes(x = sal_pss, y = ta_umolkg)) +
  geom_point()

# Goleta - low flow
# NO distinct pattern
full_fresh %>% 
  filter(Site == "Goleta") %>% 
  ggplot(aes(x = t_C, y = pH_total)) +
  geom_point()

# no real pattern - non-linear (evaluate freshwater signal)
full_fresh %>% 
  filter(Site == "Goleta") %>% 
  ggplot(aes(x = sal_pss, y = ta_umolkg)) +
  geom_point()

# W of Gaviota SP - very low flow
# NO distinct pattern
full_fresh %>% 
  filter(Site == "W of Gaviota SP") %>% 
  ggplot(aes(x = t_C, y = pH_total)) +
  geom_point()

# small salinity - alkalinity trend
full_fresh %>% 
  filter(Site == "W of Gaviota SP") %>% 
  ggplot(aes(x = sal_pss, y = ta_umolkg)) +
  geom_point()

# Columbia River

full_fresh %>% 
  filter(Site == "Columbia River") %>% 
  ggplot(aes(x = t_C, y = pH_total)) +
  geom_point()

# small salinity - alkalinity trend
full_fresh %>% 
  filter(Site == "Columbia River") %>% 
  ggplot(aes(x = sal_pss, y = ta_umolkg)) +
  geom_point()

# Tester

full_fresh %>% 
  filter(Site == "Ventura; Ventura River") %>% 
  ggplot(aes(x = t_C, y = pH_total)) +
  geom_point()

# small salinity - alkalinity trend
full_fresh %>% 
  filter(Site == "Ventura; Ventura River") %>% 
  ggplot(aes(x = sal_pss, y = ta_umolkg)) +
  geom_point()

# Use Tomales Bay, Refugio, and Columbia River
# Other compelling cases: Umpqua River, Neskowin Creek, Rogue River, Tillamook Bay (Heitmiller Creek; Watseco Creek), Stillwater Cove (Stockhoff Creek), San Vicente Creek
# Ventura has the opposite pattern (inverse TA/salinity), but barely

# Add SF Bay to the mix. There are a TON of SF Bay sites so need to filter by lat/lon
SF_Bay$Lat <- as.numeric(SF_Bay$Lat)
SF_Bay$Lon <- as.numeric(SF_Bay$Lon)
SF_Bay <- full2[(full2$Lat >= 37.4 & full2$Lat <= 38.2705 & full2$Lon <= -122.515255 & full2$Lon >= -121.704),]

# Try plotting temp and pH for SF Bay - there is no pH or TA data, only temp and salinity
SF_Bay %>% 
  ggplot(aes(x = t_C, y = pH_total)) +
  geom_point()

sum(is.na(SF_Bay$ta_umolkg)) # Confirmed - NA values for TA (and pH) for the entire SF Bay dataset

# Columbia River (offshore points)
CR_lat <- c(46.857,46.9859,47.1336)
CR_lon <- c(-124.244, -124.566, -124.272)

CR <- full2 %>% filter(Lat %in% CR_lat) %>% filter(Lon %in% CR_lon)

sum(!is.na(CR$pH_total)) # pH data is limited
sum(!is.na(CR$t_C)) # temp data is more common

View(CR[(!is.na(CR$pH_total)),]) # However times with pH data do not have temp (parameters do not overlap)

CR_coast <- full_fresh %>% filter(Site == "Columbia River")



