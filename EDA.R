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
library(lme4)

### Sorting datasets for examining freshwater impact on OA in California, Washington, Oregon.
# Filter by datasets with time series inside relevant dates
# full is loaded from RData `clean_all_distance_offshore_seacarb`
# Join datasets
full <- clean_all

Columbia1 <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/66_Final_OOI_SLHNov2022.csv")
OOI <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/66_Final_OOI_SLHNov2022.csv") %>% dplyr::select(-X)# OOI data (dataset 66)
Trinidad <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/67_TrinidadHeadLine/67_TrinidadHead_Nov2022SLH.csv") %>% dplyr::select(-X)
Newport <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/68_NewportHydrographicLine/68_NewportHydrographicLine_Final.csv") %>% dplyr::select(-c(X,X.1))
Columbia <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/70_CMOPSaturn2/70_CMOPSaturn2_final.csv") %>% dplyr::select(-X)
Scripps <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/71_ScrippsLaJolla/71_ScrippsLaJolla_final.csv") %>% dplyr::select(-X)
variable.names(full)
variable.names(Scripps)

# 5, 23, 21

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
coords_full <- full2 %>% filter(!is.na(ta_umolkg)) %>% distinct(latitude, longitude, dataset_id, min_time, max_time)
# Filter for datasets that are over 3 months (90 days) in time range
coords_full$time_diff <- coords_full$max_time - coords_full$min_time # Difference between youngest and oldest time
# Remove CalCOFI (dataset 25) because not nearshore
coords_ts <- coords_full %>% filter(dataset_id != 25) %>% na.omit()
# Clean up errors in coordinates. There are some crazy values with inconsistent spatial extents that alter final mapping
coords_ts <- coords_ts %>% filter(latitude < 50)
  
mcoords_ts <- st_as_sf(coords_ts, coords = c("longitude", "latitude"), crs = 4326)

st_write(mcoords_ts, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/all_sites.2March2023.kml", driver = "kml")
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
full_fresh <- full2 %>% full_join(fresh, by = c("Lat" = "Lat", "Lon" = "Lon", "Dataset_ID" = "Dataset_ID"), multiple = "all")
full_fresh$Distance <- as.numeric(full_fresh$Distance)
# Join additional Columbia sites from WCOA cruises. Assign "Columbia River" as site name so that these cruise data are included.
# full_fresh$Site <- ifelse((full_fresh$Lat >= 46) & (full_fresh$Lat <= 46.3) & (full_fresh$Lon >= -124.22) & (full_fresh$depth_m <= 10) & (!is.na(full_fresh$pH_total)), "Columbia River", full_fresh$Site)
# full_fresh$distance_offshore <- ifelse((full_fresh$Lat >= 46) & (full_fresh$Lat <= 46.3) & (full_fresh$Lon >= -124.22) & (full_fresh$Dataset_ID == 21) & (!is.na(full_fresh$pH_total)), "", full_fresh$distance_offshore)

# Filter by sites with complete values for temp and pH to check for upwelling signal (or departures thereof)
full_fresh_pH_temp <- full_fresh_pH %>% 
  filter(!is.na(t_C))

# Tally number of observations per Site with temp and pH data
View(full_fresh_pH_temp %>% group_by(Site) %>% tally())

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
  filter(Site == "Vandenberg State Marine Reserve (Santa Ynez River)") %>% 
  ggplot(aes(x = t_C, y = pH_total)) +
  geom_point()

# Refugio - very low flow
# SMALL pattern from upwelling
full_fresh %>% 
  filter(Site == "Vandenberg State Marine Reserve (Santa Ynez River)") %>% 
  ggplot(aes(x = sal_pss, y = ta_umolkg)) +
  geom_point()

# Moderate freshwater - alkalinity pattern
full_fresh %>% 
  filter(Site == "Refugio") %>% 
  ggplot(aes(x = t_C, y = pH_total)) +
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
  filter(Site == "Eel River and Humboldt Bay") %>% 
  ggplot(aes(x = sal_pss, y = ta_umolkg)) +
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
  ggplot(aes(x = t_C, y = pH_total)) +
  geom_point()

# Tester

full_fresh %>% 
  filter(Site == "Santa Barbara Pier; Mission Creek and Sycamore Creek") %>% 
  ggplot(aes(x = t_C, y = pH_total)) +
  geom_point()

# small salinity - alkalinity trend
full_fresh %>% 
  filter(Site == "Rogue River" | Site == "Umpqua River" | Site == "Neskowin Creek" | Site == "Tillamook Bay; Heitmiller Creek; Watseco Creek") %>% 
  ggplot(aes(x = sal_pss, y = ta_umolkg)) +
  geom_point()

# Use Tomales Bay, Refugio, and Columbia River, collection of SoCal
# Other compelling cases: Umpqua River, Neskowin Creek, Rogue River, Tillamook Bay; Heitmiller Creek; Watseco Creek, Stillwater Cove (Stockhoff Creek), San Vicente Creek
# Ventura has the opposite pattern (inverse TA/salinity), but barely

###### Add stream gauge data where possible

full_fresh$date <- as.Date(full_fresh$time_utc)

cond_fresh <- full_fresh %>% filter(!is.na(ta_umolkg)) %>% filter(!is.na(Site)) %>% 
  group_by(date, Site, Dataset_ID) %>% summarize(pH_daily = mean(pH_total, na.rm = TRUE),
                                              temp_daily = mean(t_C, na.rm = TRUE),
                                              alk_daily = mean(ta_umolkg, na.rm = TRUE),
                                              sal_daily = mean(sal_pss, na.rm = TRUE),
                                              depth_mean = mean(depth_m, na.rm = TRUE),
                                              distance = mean(Distance, na.rm = TRUE))

### Tomales Bay
# Tomales Bay - determine relevant dates
Tomales <- cond_fresh %>% filter(Site == "Tomales Bay")
# Min time = 2012-08-21 18:17:00
# Max time = 2019-07-25 22:24:00

# Match stream gauge data by date
Tomales_stream1 <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/Data/flow/TomalesBay_Lagunitas_USGS.csv")
Tomales_stream1$date <- mdy(Tomales_stream1$datetime)
Tomales_cond <- Tomales %>% left_join(Tomales_stream1, by = "date")

Tomales_stream2 <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/Data/flow/TomalesBay_Walker_USGS.csv")
Tomales_stream2$date <- mdy(Tomales_stream2$datetime)
Tomales_cond <- Tomales_cond %>% left_join(Tomales_stream2, by = "date")
  
Tomales_cond$discharge_cfs <- Tomales_cond$discharge_cfs.x  + Tomales_cond$discharge_cfs.y

Tomales_cond %>% ggplot(aes(x = log(discharge_cfs), y = alk_daily)) + geom_point() + theme_bw() + geom_smooth(method = "lm")
Tomales_cond %>% ggplot(aes(x = log(discharge_cfs), y = pH_daily)) + geom_point() + theme_bw()


### Columbia River
# Filter out dataset 70 (no pH data)
ColumbiaR <- cond_fresh %>% filter(Site == "Columbia River") %>% filter(Dataset_ID != 70)
# Min time = 2010-09-07 07:00:00
# Max time = 2016-05-30 07:00:00
# Match Columbia River dataset to stream gauge data

Columbia_stream <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/Data/flow/ColumbiaRiver_USGS.csv")
Columbia_stream$date <- mdy(Columbia_stream$datetime)
Columbia_cond <- ColumbiaR %>% left_join(Columbia_stream, by = "date")

Columbia_cond %>% ggplot(aes(x = log(discharge_tidalf_cfs), y = alk_daily)) + geom_point() + theme_bw() + geom_smooth(method = "lm")

### Rogue, Umpqua, Neskowin, Tillamook, Chetco
Rogue <- cond_fresh %>% filter(Site == "Rogue River")
# Min time = 2010-09-02 07:00:00
# Max time = 2015-05-08 07:00:00

# Match OA dataset to stream gauge data
Rogue_stream <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/Data/flow/Rogue_USGS.csv")
Rogue_stream$date <- mdy(Rogue_stream$datetime)
Rogue_cond <- Rogue %>% left_join(Rogue_stream, by = "date")

Rogue_cond %>% ggplot(aes(x = discharge_cfs, y = alk_daily)) + geom_point() + geom_smooth(method = "lm")

### Rogue, Umpqua, Neskowin, Tillamook, Chetco
Umpqua <- cond_fresh %>% filter(Site == "Umpqua River")
# Min time = 2010-09-02 07:00:00
# Max time = 2015-05-08 07:00:00

# Match OA dataset to stream gauge data
Umpqua_stream <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/Data/flow/Umpqua_USGS.csv")
Umpqua_stream$date <- mdy(Umpqua_stream$datetime)
Umpqua_cond <- Umpqua %>% left_join(Umpqua_stream, by = "date")

Umpqua_cond %>% ggplot(aes(x = discharge_cfs, y = alk_daily)) + geom_point() + geom_smooth(method = "lm")

### Rogue, Umpqua, Neskowin, Tillamook
Tillamook <- cond_fresh %>% filter(Site == "Tillamook Bay; Heitmiller Creek; Watseco Creek")
# Min time = 2010-09-02 07:00:00
# Max time = 2015-05-08 07:00:00

# Match OA dataset to stream gauge data
Tillamook_stream1 <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/Data/flow/Tillamook.Trask_USGS.csv")
Tillamook_stream1$date <- mdy(Tillamook_stream1$datetime)
Tillamook_stream1 <- Tillamook_stream1 %>% dplyr::select(c(site_no, date, discharge_cfs))
Tillamook_cond <- Tillamook %>% left_join(Tillamook_stream1, by = "date")

Tillamook_stream2 <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/Data/flow/Tillamook.Wilson_USGS.csv")
Tillamook_stream2$date <- mdy(Tillamook_stream2$datetime)
Tillamook_stream2 <- Tillamook_stream2 %>% dplyr::select(c(site_no, date, discharge_cfs))
Tillamook_cond <- Tillamook_cond %>% left_join(Tillamook_stream2, by = "date")

Tillamook_cond$discharge_cfs <- Tillamook_cond$discharge_cfs.x  + Tillamook_cond$discharge_cfs.y

Tillamook_cond %>% ggplot(aes(x = discharge_cfs, y = alk_daily)) + geom_point() + geom_smooth(method = "lm")

### Rogue, Umpqua, Neskowin, Tillamook, Chetco
Chetco <- cond_fresh %>% filter(Site == "Chetco River")
# Join stream data
Chetco_stream <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/Data/flow/Chetco_USGS.csv")
Chetco_stream$date <- mdy(Chetco_stream$datetime)
Chetco_stream <- Chetco_stream %>% dplyr::select(c(site_no, date, discharge_cfs))
Chetco_cond <- Chetco %>% left_join(Chetco_stream, by = "date")

Chetco_cond %>% ggplot(aes(x = discharge_cfs, y = alk_daily)) + geom_point() + geom_smooth(method = "lm")


### Combine all rivers above to view pattern
Till_bind <- Tillamook_cond %>% dplyr::select(c(date, pH_daily, temp_daily, alk_daily, sal_daily, discharge_cfs, distance))
Till_bind$Site <- "Tillamook"
Rogue_bind <- Rogue_cond %>% dplyr::select(c(date, pH_daily, temp_daily, alk_daily, sal_daily, discharge_cfs, distance))
Rogue_bind$Site <- "Rogue"
Umpqua_bind <- Umpqua_cond %>% dplyr::select(c(date, pH_daily, temp_daily, alk_daily, sal_daily, discharge_cfs, distance))
Umpqua_bind$Site <- "Umpqua"
Chetco_bind <- Chetco_cond %>% dplyr::select(c(date, pH_daily, temp_daily, alk_daily, sal_daily, discharge_cfs, distance))
Chetco_bind$Site <- "Chetco"

NW <- do.call("rbind", list(Till_bind, Rogue_bind, Umpqua_bind, Chetco_bind))
NW %>% ggplot(aes(x = log(discharge_cfs), y = alk_daily, group = Site)) + geom_point(aes(col = Site)) + geom_smooth(aes(color = Site), method = "lm") + facet_grid(rows = vars(Site)) + theme_bw() + ylim(1400, 2375)
NW %>% ggplot(aes(x = log(discharge_cfs), y = pH_daily, group = Site)) + geom_point(aes(col = Site)) + facet_grid(rows = vars(Site)) + theme_bw()

Columbia_cond %>% ggplot(aes(x = log(discharge_cfs), y = alk_daily)) + geom_point() + geom_smooth(method = "lm")
Columbia_cond %>% ggplot(aes(x = log(discharge_cfs), y = pH_daily)) + geom_point() + geom_smooth(method = "lm") + theme_bw()

# 0 - 10; 10 - 50

### SoCal

# Goleta
# Match stream gauge data by date
Goleta <- cond_fresh %>% filter(Site == "Goleta")
# Join stream data
Goleta_stream <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/Data/flow/Goleta_USGS.csv")
Goleta_stream$date <- mdy(Goleta_stream$datetime)
Goleta_stream <- Goleta_stream %>% dplyr::select(c(site_no, date, discharge_cfs))
Goleta_cond <- Goleta %>% left_join(Goleta_stream, by = "date")

Goleta_cond %>% ggplot(aes(x = log(discharge_cfs), y = alk_daily)) + geom_point()
Goleta_cond %>% ggplot(aes(x = log(discharge_cfs), y = pH_daily)) + geom_point()

# Ventura
Ventura <- cond_fresh %>% filter(Site == "Ventura; Ventura River")
# min date = 2010-09-07
# max date = 2015-05-06
Ventura_stream <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/Data/flow/Ventura_USGS.csv")
Ventura_stream$date <- mdy(Ventura_stream$datetime)
Ventura_stream <- Ventura_stream %>% dplyr::select(c(site_no, date, discharge_cfs))
Ventura_cond <- Ventura %>% left_join(Ventura_stream, by = "date")

Ventura_cond %>% ggplot(aes(x = log(discharge_cfs), y = alk_daily)) + geom_point()
Ventura_cond %>% ggplot(aes(x = log(discharge_cfs), y = pH_daily)) + geom_point()


# Santa Barbara Pier; Mission Creek and Sycamore Creek
SB <- cond_fresh %>% filter(Site == "Santa Barbara Pier; Mission Creek and Sycamore Creek")

SB_stream <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/Data/flow/SantaBarbara_USGS.csv")
SB_stream$date <- mdy(SB_stream$datetime)
SB_stream <- SB_stream %>% dplyr::select(c(site_no, date, discharge_cfs))
SB_cond <- SB %>% left_join(SB_stream, by = "date")

SB_cond %>% ggplot(aes(x = discharge_cfs, y = alk_daily)) + geom_point() + theme_bw()
SB_cond %>% ggplot(aes(x = log(discharge_cfs), y = pH_daily)) + geom_point()


##### GLMM
Columbia_cond2 <- Columbia_cond %>% dplyr::select(c(Site, date, pH_daily, temp_daily, alk_daily, sal_daily, discharge_cfs, distance)) %>% ungroup() %>% sample_n(10)
Tomales_cond2 <- Tomales_cond %>% dplyr::select(c(Site, date, pH_daily, temp_daily, alk_daily, sal_daily, discharge_cfs, distance)) %>% ungroup() %>% sample_n(10)
Goleta_cond2 <- Goleta_cond %>% dplyr::select(c(Site, date, pH_daily, temp_daily, alk_daily, sal_daily, discharge_cfs, distance)) %>% ungroup() %>% sample_n(10)
Ventura_cond2 <- Ventura_cond %>% dplyr::select(c(Site, date, pH_daily, temp_daily, alk_daily, sal_daily, discharge_cfs, distance))
SB_cond2 <- SB_cond %>% dplyr::select(c(Site, date, pH_daily, temp_daily, alk_daily, sal_daily, discharge_cfs, distance)) %>% ungroup() %>% sample_n(10)

cond2 <- do.call("rbind", list(Columbia_cond2, NW, Tomales_cond2, Goleta_cond2, SB_cond2))

cond2$year <- year(cond2$date)
cond2$discharge_scale <- scale(cond2$discharge_cfs)
cond2$distance_scale <- scale(cond2$distance)

m <- lmer(alk_daily ~ discharge_scale + distance_scale + year + ( 1| Site), data = cond2)
m2 <- lm(alk_daily ~ discharge_scale + distance_scale + year + Site, data = cond2)

summary(m)
summary(m2)


### List to search inland
# Goleta is ideal
# Cambria/Hearst State Park (Santa Rosa Creek and Leffingwell Creek; Crystal Cove (2); San Clamente

###### GLORICH dataset
glo <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/GLORICH/hydrochemistry.csv")
glo <- glo %>% filter(!is.na(Alkalinity))
coords <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/GLORICH/sampling_locations.csv")
glo <- left_join(glo, coords, by = "STAT_ID")
glo <- glo %>% filter(Country == "USA") %>% filter(State == "OR" | State == "WA" | State == "CA" | State == "OREGON" | State == "WASHINGTON" | State == "CALIFORNIA")
hist(glo$Alkalinity)





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



