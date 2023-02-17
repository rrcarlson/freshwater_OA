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

# Bind datasets
OOI$distance_offshore <- NA # Create a blank column for distance_offshore as this is missing in individual datasets
Trinidad$distance_offshore <- NA
Newport$distance_offshore <- NA
Columbia$distance_offshore <- NA
Scripps$distance_offshore <- NA

full2 <- do.call("rbind", list(full, OOI, Trinidad, Newport, Columbia, Scripps)) # Bind datasets

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

# The above seems to be missing a lot - Recalculate time ranges and and plot full data
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
full2$latitude <- sub("0+$", "", as.character(full2$latitude))
full2$longitude <- sub("0+$", "", as.character(full2$longitude))

# Import sites tagged for river proximity
fresh <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/sites.csv")
fresh$Lat <- sub("0+$", "", as.character(fresh$Lat))
fresh$Lon <- sub("0+$", "", as.character(fresh$Lon))
fresh <- fresh %>% dplyr::select(c("Site", "Dataset_ID","Distance","Reach_ID","Flow_range_midFeb","Mean_flow","Lat","Lon"))
full2 <- full2 %>% rename("Lat" = "latitude", "Lon" = "longitude", "Dataset_ID" = "dataset_id") # There were 70ish mismatched dataset IDs in first try at joining these datasets. Ensure datest ID in tagged freshwater dataset is correct.

# Join reach, site name, flow, and distance information to relevant sites
full_fresh <- full2 %>% right_join(fresh, by = c("Lat" = "Lat", "Lon" = "Lon", "Dataset_ID" = "Dataset_ID"), multiple = "all")
View(full_fresh %>% filter(Site == "Eel River and Humboldt Bay"))
full_fresh_pH <- full_fresh %>% filter(!is.na(pH_total))

