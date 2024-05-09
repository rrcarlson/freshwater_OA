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
library(foreign) # allows you to read DBF
library(rgdal)
library(readr)

### Sorting datasets for examining freshwater impact on OA in California, Washington, Oregon.
# Filter by datasets with time series inside relevant dates
# full is loaded from RData `clean_all_distance_offshore_seacarb`
# Join datasets
full <- clean_all
variable.names(full)

# Calculate time ranges
full <- full %>% group_by(dataset_id) %>% mutate(min_time = min(time_utc, na.rm = TRUE),
                                                   max_time = max(time_utc, na.rm = TRUE))
# Identify unique coordinates (unique points in dataset) for mapping purposes
coords_full <- full %>% filter(!is.na(ta_umolkg)) %>% distinct(latitude, longitude, dataset_id, min_time, max_time)
# Filter for datasets that are over 3 months (90 days) in time range
coords_full$time_diff <- coords_full$max_time - coords_full$min_time # Difference between youngest and oldest time
# Remove CalCOFI (dataset 25) because not nearshore
coords_ts <- coords_full %>% filter(dataset_id != 25) %>% na.omit()
# Clean up errors in coordinates. There are some crazy values with inconsistent spatial extents that alter final mapping
coords_ts <- coords_ts %>% filter(latitude < 50)
  
mcoords_ts <- st_as_sf(coords_ts, coords = c("longitude", "latitude"), crs = 4326)

st_write(mcoords_ts, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/all_sites.2March2023.kml", driver = "kml")
st_write(mcoords_ts,"/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/all_sites.9Feb2023.shp")

# Extract only points where alkalinity is nonzero.
full$latitude <- as.numeric(sub("0+$", "", as.character(full$latitude))) %>% round(digits = 4)
full$longitude <- as.numeric(sub("0+$", "", as.character(full$longitude))) %>% round(digits = 4)

full <- full %>% filter(!is.na(ta_umolkg)) # Reduces dataset from ~14 million to 12152 sites
full$date <- date(full$time_utc)
full_trim <- full %>% dplyr::select(c(dataset_id, latitude, longitude, depth_m, time_utc, date, t_C, sal_pss, pH_total, tCO2_umolkg, ta_umolkg, ta_flag, do_umolkg, chl_ugL, si_umolkg, nh4_umolkg, no3_umolkg, no2_umolkg, po3_umolkg, min_time, max_time, distance_offshore))
# Several observations appear to have repeated observations per day over several hours. Find mean per day.
full_daily <- full_trim %>% group_by(date, latitude, longitude, depth_m, dataset_id) %>% summarize(t_C = mean(t_C, na.rm = TRUE),
                                                                                              sal_pss = mean(sal_pss, na.rm = TRUE),
                                                                                              pH_total = mean(pH_total, na.rm = TRUE),
                                                                                              tCO2_umolkg = mean(tCO2_umolkg, na.rm = TRUE),
                                                                                              ta_umolkg = mean(ta_umolkg, na.rm = TRUE),
                                                                                              do_umolkg = mean(do_umolkg, na.rm = TRUE),
                                                                                              chl_ugL = mean(chl_ugL, na.rm = TRUE),
                                                                                              si_umolkg = mean(si_umolkg, na.rm = TRUE),
                                                                                              nh4_umolkg = mean(nh4_umolkg, na.rm = TRUE),
                                                                                              no3_umolkg = mean(no3_umolkg, na.rm = TRUE),
                                                                                              no2_umolkg = mean(no2_umolkg, na.rm = TRUE),
                                                                                              po3_umolkg = mean(po3_umolkg, na.rm = TRUE),
                                                                                              distance_offshore = mean(distance_offshore, na.rm = TRUE))
# full_daily_ca <- full_daily %>% dplyr::filter(latitude > 32.48 & latitude < 42.02) %>% st_as_sf(coords = c("longitude","latitude"), crs = 4326)
# full_daily_or <- full_daily %>% dplyr::filter(latitude > 41.998 & latitude < 46.08) %>% st_as_sf(coords = c("longitude","latitude"), crs = 4326)
# full_daily_wa <- full_daily %>% dplyr::filter(latitude > 46.08 & latitude < 49.174) %>% st_as_sf(coords = c("longitude","latitude"), crs = 4326)
# full_daily_pnw <- full_daily %>% dplyr::filter(latitude > 41.998 & latitude < 49.174) %>% st_as_sf(coords = c("longitude","latitude"), crs = 4326)
full_daily <- full_daily %>% st_as_sf(coords = c("longitude","latitude"), crs = 4326) # Convert to shapefile


########## Inland alkalinity mapping

###### GLORICH dataset
glo <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Alkalinity/GLORICH/hydrochemistry.csv")
glo <- glo %>% filter(!is.na(Alkalinity)) # Filter GLORICH for only rows where alkalinity data exists
# Join GLORICH hydrochem data (i.e., alkalinity) to GLORICH station log (contains lat/lon)
coords <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Alkalinity/GLORICH/sampling_locations.csv")
glo <- left_join(glo, coords, by = "STAT_ID")
# Filter to only stations inside US and West Coast
glo <- glo %>% filter(Country == "USA") %>% filter(State == "OR" | State == "WA" | State == "CA" | State == "OREGON" | State == "WASHINGTON" | State == "CALIFORNIA")
hist(glo$Alkalinity)
# Write out spatial dataset
write.csv(glo, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Alkalinity/GLORICH/hydrochem_joined.csv")

###### SWAMP (State of CA) dataset
swamp <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Alkalinity/SWAMP/SWAMP_original.csv")
swamp$ta_umolkg <- swamp$Result * 10 * 2
hist(swamp$ta_umolkg)
swamp_sf <- st_as_sf(swamp, coords = c("TargetLongitude", "TargetLatitude"))
st_write(swamp_sf, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Alkalinity/SWAMP/SWAMP_original.shp")

# The following was cleaned in Excel prior to use. Raw data is messy in an inconsistent way - hard to clean/reshape in R.
new_usgs <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Alkalinity/USGS/USGS_cleaned.csv", header = TRUE)
new_usgs <- new_usgs %>% spread(parm_cd,result_va)
colnames(new_usgs)[11:17] <- c("c00418","c29801","c29802","c29803","c39036","c39086","c39087")

new_usgs %>% group_by(site_no) %>% tally() %>% nrow() # 6601 sites in California, Washington, and Oregon

new_usgs <- new_usgs %>% mutate(alkalinity_mgL = rowMeans(new_usgs[,11:17], na.rm = TRUE))
new_usgs <- new_usgs %>% mutate(alkalinity_umolkg = 2*10*alkalinity_mgL)

# Ensure that all rows have an alkalinity value; all rows in this dataset do so the USGS data retrieval (for nonzero alkalinity) was correct
new_usgs <- new_usgs %>% filter(!is.na(alkalinity_mgL))

# Remove one anomalously high value
new_usgs <- new_usgs[-16380,]

# Load site information to join USGS alkalinity with lat/lon data
new_usgs_sites <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Alkalinity/USGS/USGS_sites.csv") %>% dplyr::select(c(site_no, site_nom, Lat, Lon))
new_usgs <- left_join(new_usgs, new_usgs_sites, by = "site_no")
new_usgs$date <- as.Date(new_usgs$sample_dt, "%m/%d/%y")
new_usgs <- new_usgs %>% dplyr::select(-sample_dt)

write.csv(new_usgs, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Alkalinity/USGS/USGS_joined.csv")

# Finish reformatting USGS and SWAMP datasets so that they can be bounded
usgs_slim <- new_usgs %>% dplyr::select(c(site_no, alkalinity_umolkg, site_nom, Lat, Lon, date))
colnames(usgs_slim)[2] <- "ta_umolkg"
swamp_slim <- swamp %>% dplyr::filter(DataQuality != "Extensive review needed") %>% dplyr::select(c(StationCode, ta_umolkg, StationName, TargetLatitude, TargetLongitude, SampleDate, DataQuality))
swamp_slim$date <- as.Date(swamp_slim$SampleDate)
swamp_slim <- swamp_slim %>% dplyr::select(-c(SampleDate, DataQuality)) %>% dplyr::filter(!is.na(ta_umolkg))
colnames(swamp_slim) <- colnames(usgs_slim)

# This is all inland alkalinity data from USGS 2007 - 2023 and SWAMP (California state)
ta_land <- rbind(swamp_slim, usgs_slim)
ta_land <- ta_land %>% filter(!is.na(Lon))

####### Aggregate data by watershed (HUC10 and HUC12)
# Create single layer of all HUC12s
huc12_ca <- st_read("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/CA/HUC12/ACE_HUC12s_WebMerc_1mXY.shp") %>% dplyr::select(HUC12, geometry)
huc12_or <- st_read("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/OR/HUC12/hydrologic_units_wbdhu12_a_or.shp") %>% dplyr::select(huc12, geometry)
huc12_wa <- st_read("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/WA/HUC12/hydrologic_units_wbdhu12_a_wa.shp") %>% dplyr::select(huc12, geometry)

colnames(huc12_or) <- c("HUC12", "geometry")
colnames(huc12_wa) <- c("HUC12", "geometry")

huc12 <- rbind(huc12_ca, huc12_or, huc12_wa)
st_write(huc12, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/All/HUC12/huc12_Wcoast.shp")

# Create single layer of all HUC10s
huc10_ca <- st_read("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/CA/HUC10/WBD_USGS_HUC10_CA.shp", crs = 6414) %>% dplyr::select(HUC10, geometry) %>% st_transform(crs = 4326)
huc10_or <- st_read("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/OR/HUC10/hydrologic_units_wbdhu10_a_or.shp") %>% dplyr::select(huc10, geometry)
huc10_wa <- st_read("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/WA/HUC10/hydrologic_units_wbdhu10_a_wa.shp") %>% dplyr::select(huc10, geometry)

colnames(huc10_or) <- c("HUC10", "geometry")
colnames(huc10_wa) <- c("HUC10", "geometry")

huc10 <- rbind(huc10_ca, huc10_or, huc10_wa)
st_write(huc10, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/All/HUC10/huc10_Wcoast.shp")

# Create single layer of all HUC8s
huc8_ca <- st_read("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/CA/HUC8/hydrologic_units_wbdhu8_a_ca.shp") %>% dplyr::select(huc8, geometry)
huc8_or <- st_read("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/OR/HUC8/hydrologic_units_wbdhu8_a_or.shp") %>% dplyr::select(huc8, geometry)
huc8_wa <- st_read("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/WA/HUC8/hydrologic_units_wbdhu8_a_wa.shp") %>% dplyr::select(huc8, geometry)

colnames(huc8_ca) <- c("HUC8", "geometry")
colnames(huc8_or) <- c("HUC8", "geometry")
colnames(huc8_wa) <- c("HUC8", "geometry")

huc8 <- rbind(huc8_ca, huc8_or, huc8_wa)
st_write(huc8, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/All/HUC8/huc8_Wcoast.shp")

# Convert USGS + SWAMP alkalinity points to sf object
ta_land_shp <- st_as_sf(ta_land, coords = c("Lon","Lat"), crs = 4326)

# Join points to watershed
sf_use_s2(FALSE)
ta_in_huc8 <- st_join(ta_land_shp, huc8, join = st_within)
ta_in_huc10 <- st_join(ta_land_shp, huc10, join = st_within)
ta_in_huc12 <- st_join(ta_land_shp, huc12, join = st_within)

# Some rows are repeated for some reason, so ensure that only one row per unique values of all columns
ta_in_huc8 <- ta_in_huc8 %>% distinct()
ta_in_huc10 <- ta_in_huc10 %>% distinct()
ta_in_huc12 <- ta_in_huc12 %>% distinct()

# Mean, maximum, and mode value per HUC10 and HUC12. DO NOT DO NA.RM = TRUE because there are alkalinity values at all points and don't want to conflate 0 alkalinity with 0 points within HUC
huc8_ta <- ta_in_huc8 %>% as.data.frame() %>% group_by(HUC8) %>% summarize(mean_ta_huc8 = mean(ta_umolkg),
                                                                           max_ta_huc8 = max(ta_umolkg),
                                                                           min_ta_huc8 = min(ta_umolkg)) %>% right_join(huc8, by = "HUC8", multiple = "all") %>% dplyr::filter(!is.na(mean_ta_huc8)) %>% st_as_sf()
huc10_ta <- ta_in_huc10 %>% as.data.frame() %>% group_by(HUC10) %>% summarize(mean_ta_huc10 = mean(ta_umolkg),
                                                                              max_ta_huc10 = max(ta_umolkg),
                                                                              min_ta_huc10 = min(ta_umolkg)) %>% right_join(huc10, by = "HUC10", multiple = "all") %>% dplyr::filter(!is.na(mean_ta_huc10)) %>% st_as_sf()
huc12_ta <- ta_in_huc12 %>% as.data.frame() %>% group_by(HUC12) %>% summarize(mean_ta_huc12 = mean(ta_umolkg),
                                                                              max_ta_huc12 = max(ta_umolkg),
                                                                              min_ta_huc12 = min(ta_umolkg)) %>% right_join(huc12, by = "HUC12", multiple = "all") %>% dplyr::filter(!is.na(mean_ta_huc12)) %>% st_as_sf()
st_write(huc8_ta, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/All/HUC8/huc8_Wcoast_w_alk.shp")
st_write(huc10_ta, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/All/HUC10/huc10_Wcoast_w_alk.shp")
st_write(huc12_ta, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/WBD/All/HUC12/huc12_Wcoast_w_alk.shp")

######## Calculate, for OA dataset, watershed-based alkalinity at HUC8, HUC10, and HUC12 levels
full_daily_coast <- full_daily %>% dplyr::filter(distance_offshore <=10) # Filter OA daily dataset to those points within 10 km of coast
# Find nearest HUC8, 10, and 12 to each OA point
nearhuc8 <- st_nearest_feature(full_daily_coast, huc8_ta)
nearhuc10 <- st_nearest_feature(full_daily_coast, huc10_ta)
nearhuc12 <- st_nearest_feature(full_daily_coast, huc12_ta)
OA_join <- full_daily_coast
# We don't need to know the distance to the watershed/HUC because we have distance to nearest stream
OA_join <- cbind(OA_join, (huc8_ta[nearhuc8,] %>% as.data.frame() %>% dplyr::select(-geometry)))
OA_join <- cbind(OA_join, (huc10_ta[nearhuc10,] %>% as.data.frame() %>% dplyr::select(-geometry)))
OA_join <- cbind(OA_join, (huc12_ta[nearhuc12,] %>% as.data.frame() %>% dplyr::select(-geometry)))

st_write(OA_join, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/OA/OA_coastal_alk.shp")
OA_join <- st_read("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/OA/OA_coastal_alk.shp")

# Load NHD per region
PNW <- st_read("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/PNW/NHDPlus/Hydrography/NHDFlowline.shp")
CA <- st_read("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/CA/NHDPlus/NHDSnapshot/Hydrography/NHDFlowline.shp")
# Write a function for ingesting mean monthly flow data from individual DBFs into a dataframe (one column per month)
# Attributes of interest are as follows. We will ingest Q0001E from each monthly DBF. 'flow_pick' below can be rewritten for another attribute if desired.
#     Q0001A = Initial estimated flow at downstream end of NHD Flowline
#     Q0001B = Flow with anticipated evapotranspiration subtracted - this may be important in SoCal
#     Q0001C = Flow with Flow with Reference Gage Regression applied to Q0001B (cfs) = best estimate of natural flow according to https://www.youtube.com/watch?v=9ZJIpRPYLi4
#     Q0001D = PlusFlowAR (cfs) = Flow from Q0001C with manmade obstructions like dams or water diversions taken into account, largely crowdsourced
#     Q0001E = Flow from gage adjustment (cfs) = Best NHDPlus flow estimates for actual flows. Adjusts flows downstream of NWIS gages (which meet certain criteria) based on gage data + incremental flow estimates. NHDPlus HR network features that are upstream from gages are adjusted for observed gage-based flow. Affects estiamtes downstream from gages to better feclect flow alterations not taken into account in first four steps.

# quantile(PNW_flow$Q0001A) # Visualize the kinds of flow values we're dealing with

# Function to extract only column Q0001E from each month's DBF
flow_pick <- function(fp){
  PNW_flow <- read.dbf(fp) %>% dplyr::select(c(ComID,Q0001E))
}

# Combine Q0001E values from each month into one dataframe
# For PNW
fp_ls_pnw <- list.files("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/PNW/EROMExtension/EROM_monthly/", full.names = TRUE)
flow_ls_pnw <- lapply(fp_ls_pnw, flow_pick)
flow_permonth_pnw <- do.call(cbind.data.frame, flow_ls_pnw)
flow_permonth_pnw <- flow_permonth_pnw[,-c(3,5,7,9,11,13,15,17,19,21,23)] # Take out repetitions of "ComID"
colnames(flow_permonth_pnw) <- c("ComID","Jan","Feb","March","April","May","June","July","Aug","Sep","Oct","Nov","Dec")

# For CA
fp_ls_ca <- list.files("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/CA/EROMExtension/EROM_Monthly/", full.names = TRUE)
flow_ls_ca <- lapply(fp_ls_ca, flow_pick)
flow_permonth_ca <- do.call(cbind.data.frame, flow_ls_ca)
flow_permonth_ca <- flow_permonth_ca[,-c(3,5,7,9,11,13,15,17,19,21,23)] # Take out repetitions of "ComID"
colnames(flow_permonth_ca) <- c("ComID","Jan","Feb","March","April","May","June","July","Aug","Sep","Oct","Nov","Dec")

# Attach mean annual flow to mean monthly flow table
ma_pnw <- read.dbf("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/PNW/EROMExtension/EROM_MA0001.DBF")
flow_pnw <- ma_pnw %>% dplyr::select(c(ComID, Q0001E)) %>% full_join(flow_permonth_pnw, by = "ComID")
colnames(flow_pnw)[2] <- "ma"

ma_ca <- read.dbf("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/CA/EROMExtension/EROM_MA0001.DBF")
flow_ca <- ma_ca %>% dplyr::select(c(ComID, Q0001E)) %>% full_join(flow_permonth_ca, by = "ComID")
colnames(flow_ca)[2] <- "ma"

# Attribute flow to NHD network
colnames(PNW)[1] <- "ComID"
colnames(CA)[1] <- "ComID"
PNW <- PNW %>% left_join(flow_pnw, by = "ComID")
PNW_rivers <- PNW %>% dplyr::filter(FCODE == 46000 | FCODE == 46003 | FCODE == 46006 | FCODE == 46007 | FCODE == 55800) # feature codes for "stream/river"; see https://files.hawaii.gov/dbedt/op/gis/data/NHD%20Complete%20FCode%20Attribute%20Value%20List.pdf
PNW_coasts <- PNW %>% dplyr::filter(FCODE == 56600)

CA <- CA %>% left_join(flow_ca, by = "ComID")
CA_rivers <- CA %>% dplyr::filter(FCODE == 46000 | FCODE == 46003 | FCODE == 46006 | FCODE == 46007 | FCODE == 55800)
CA_coasts <- CA %>% dplyr::filter(FCODE == 56600)

st_write(PNW, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/PNW/NHDPlus/Hydrography/NHDFlowline_plus_flow.shp")
st_write(PNW_rivers, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/PNW/NHDPlus/Hydrography/NHDRivers_plus_flow.shp")
st_write(PNW_coasts, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/PNW/NHDPlus/Hydrography/NHDCoasts.shp")
st_write(CA, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/CA/NHDPlus/NHDSnapshot/Hydrography/NHDFlowline_plus_flow.shp")
st_write(CA_rivers, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/CA/NHDPlus/NHDSnapshot/Hydrography/NHDRivers_plus_flow.shp")
st_write(CA_coasts, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/CA/NHDPlus/NHDSnapshot/Hydrography/NHDCoasts.shp")

# Need to find only stream reaches that intersect the coast, so that when you buffer, you don't fetch inland stream reaches as well as coastal streams and double count
coasts <- rbind(PNW_coasts, CA_coasts)
rivers <- rbind(PNW_rivers, CA_rivers)

coastal_riv_df <- st_touches(st_zm(coasts), st_zm(rivers)) %>% as.data.frame()
coastal_riv_index <- unique(coastal_riv_df$col.id)
coastal_riv <- rivers[coastal_riv_index,]

# The above misses a few flowlines that end a hair's breadth from the coast. Manually identified those Flowlines, must add them back.
ls_sub <- c(23887114, 23920456, 23921184, 8317133, 3880318, 5329941, 948050160, 5329591, 17688275, 17682150, 17682488, 17682498, 17684112, 17684082, 17620018, 163864377, 17595353, 17596221, 17586992, 20365167, 20351631, 20331956)
rivers_subset <- rivers %>% dplyr::filter(ComID %in% ls_sub)
coastal_riv <- rbind(coastal_riv, rivers_subset)

st_write(coastal_riv, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/NHDcoastal_rivers.shp")
st_write(coasts, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/NHDCoastline.shp")
st_write(rivers, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/NHDRivers.shp")


# Join flow data to main frame
OA_join <- st_transform(OA_join, crs = st_crs(coastal_riv))
D <- st_distance(OA_join, coastal_riv)
nearest <- max.col(-D)

OA_join$ComID <- NA
coastal_riv_df <- coastal_riv %>% as.data.frame() %>% dplyr::select(c(ComID, FDATE, REACHCODE, FCODE, ma))
for (i in 1:length(nearest)) {
  OA_join[i,30] = coastal_riv_df[nearest[i],1]
}
OA_join <- left_join(OA_join, coastal_riv_df, by = "ComID")


# Join distance from point to the nearest river with the flow given above.
OA_join$dstnc_f_riv <- NA
for (i in 1:nrow(D)){
  OA_join$dstnc_f_riv[i] = as.numeric(D[i, nearest[i]])
}


# Find soil pH by watershed (files ingested from Google Earth Engine)
fp_ls_soilPh <- list.files("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/SoilpH/", pattern = ".csv", full.names = TRUE)
soilPh_HUC8_mean <- read.csv(fp_ls_soilPh[6]) %>% as.data.frame() %>% dplyr::select(c(HUC8, b0)) %>% rename(b0.HUC8mean = b0) %>% group_by(HUC8) %>% summarize(b0.HUC8mean = mean(b0.HUC8mean, na.rm = TRUE))
soilPh_HUC8_max <- read.csv(fp_ls_soilPh[5]) %>% as.data.frame() %>% dplyr::select(c(HUC8, b0)) %>% rename(b0.HUC8max = b0) %>% group_by(HUC8) %>% summarize(b0.HUC8max = max(b0.HUC8max, na.rm = TRUE))
soilPh_HUC10_mean <- read.csv(fp_ls_soilPh[2]) %>% as.data.frame() %>% dplyr::select(c(HUC10, b0)) %>% rename(b0.HUC10mean = b0) %>% group_by(HUC10) %>% summarize(b0.HUC10mean = mean(b0.HUC10mean, na.rm = TRUE))
soilPh_HUC10_max <- read.csv(fp_ls_soilPh[1]) %>% as.data.frame() %>% dplyr::select(c(HUC10, b0)) %>% rename(b0.HUC10max = b0) %>% group_by(HUC10) %>% summarize(b0.HUC10max = max(b0.HUC10max, na.rm = TRUE))
soilPh_HUC12_mean <- read.csv(fp_ls_soilPh[4]) %>% as.data.frame() %>% dplyr::select(c(HUC12, b0)) %>% rename(b0.HUC12mean = b0) %>% group_by(HUC12) %>% summarize(b0.HUC12mean = mean(b0.HUC12mean, na.rm = TRUE))
soilPh_HUC12_max <- read.csv(fp_ls_soilPh[3]) %>% as.data.frame() %>% dplyr::select(c(HUC12, b0)) %>% rename(b0.HUC12max = b0) %>% na.omit() %>% group_by(HUC12) %>% summarize(b0.HUC12max = max(b0.HUC12max, na.rm = TRUE))
# Join soil pH by watershed to watershed dataset
soilPh_HUC8_mean$HUC8 <- as.character(soilPh_HUC8_mean$HUC8)
soilPh_HUC8_max$HUC8 <- as.character(soilPh_HUC8_max$HUC8)
soilPh_HUC10_mean$HUC10 <- as.character(soilPh_HUC10_mean$HUC10)
soilPh_HUC10_max$HUC10 <- as.character(soilPh_HUC10_max$HUC10)
soilPh_HUC12_mean$HUC12 <- as.character(soilPh_HUC12_mean$HUC12)
soilPh_HUC12_max$HUC12 <- as.character(soilPh_HUC12_max$HUC12)

OA_join2 <- left_join(OA_join, soilPh_HUC8_mean, by = "HUC8", relationship = "many-to-many")
OA_join2 <- left_join(OA_join2, soilPh_HUC8_max, by = "HUC8", relationship = "many-to-many")
OA_join2 <- left_join(OA_join2, soilPh_HUC10_mean, by = "HUC10", relationship = "many-to-many")
OA_join2 <- left_join(OA_join2, soilPh_HUC10_max, by = "HUC10", relationship = "many-to-many")
OA_join2 <- left_join(OA_join2, soilPh_HUC12_mean, by = "HUC12", relationship = "many-to-many")
OA_join2 <- left_join(OA_join2, soilPh_HUC12_max, by = "HUC12", relationship = "many-to-many")


# Load in upwelling index
upwell <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Upwell/CUTI_daily.csv")
upwellN <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Upwell/BEUTI_daily.csv")
# Extract coordinates from point geometry in OA_join master dataset
OA_join2$lat <- st_coordinates(OA_join2)[,2]
OA_join2 <- OA_join2 %>% dplyr::filter(lat >= 32.0)
upwell <- upwell %>% gather("lat","upwell", 4:20)
upwellN <- upwellN %>% gather("lat","upwellN", 4:20)
upwell$lat_int <- parse_number(upwell$lat)
upwellN$lat_int <- parse_number(upwellN$lat)
OA_join2$lat_int <- round(OA_join2$lat) # round to nearest integer so can join to upwelling lat
OA_join2$year <- year(OA_join2$date)
OA_join2$month <- month(OA_join2$date)
OA_join2$day <- day(OA_join2$date)
upwell <- upwell %>% dplyr::select(-lat)
upwellN <- upwellN %>% dplyr::select(-lat)
# Join upwelling index to main dataframe
OA_join2 <- left_join(OA_join2, upwell, by = c("year","month","day","lat_int"))
OA_join2 <- left_join(OA_join2, upwellN, by = c("year","month","day","lat_int"))


# Extract flow rate by month, by lag month, and by season
colnames(coastal_riv)[16:27] <- c("January","February","March","April","May","June","July","August","September","October","November","December")
coastal_riv_gather <- coastal_riv %>% gather("month","month_cfs", 16:27)
coastal_riv_gather$month <- match(coastal_riv_gather$month, month.name)
coastal_riv_gather <- coastal_riv_gather %>% as.data.frame() %>% dplyr::select(c("ComID","month","month_cfs"))
OA_join2 <- left_join(OA_join2, coastal_riv_gather, by = c("month", "ComID"))


# Convert distance from shore from km to m
OA_join2$dstnc_f <- OA_join2$dstnc_f*1000
plot(OA_join2$dstnc_f, OA_join2$dstnc_f_riv) # Make sure there's a reasonable increase in distance from nearest river mouth with distance from shore. There is!


# Write dataset to file
st_write(OA_join2, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/OA/OA_master_9May2024.shp")
write.csv(OA_join2 %>% as.data.frame() %>% dplyr::select(-geometry), "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/OA/OA_master_9May2024.csv")































































####################### OLD: DO NOT USE!!!!! #######################
####################### OLD: DO NOT USE!!!!! #######################
####################### OLD: DO NOT USE!!!!! #######################
####################### OLD: DO NOT USE!!!!! #######################
####################### OLD: DO NOT USE!!!!! #######################

######## Calculate distance between each OA point and stream
# Oregon and Washington stream/river data comes from NHD. California has its own, cleaned-up version of NHD.
# As a result, Oregon and Washington have NHD layers that are disaggregated into 3-4 layers per state.
# NEED TO TRANSFORM, FILTER, AND COMBINE DISAGGREGATED NHD FILES INTO ONE COHESIVE SF OBJECT PER STATE.

# Define function for converting a list of shapefile filepaths to a list of sf objects cropped to the coast
# fpath = shapefile path; box_name is bounding box for cropping, specific to each state
convert_nhd <- function(fpath, box_name) {
  nhd_sf <- st_read(fpath) %>% 
    st_transform(crs = 4326) %>% st_zm() %>% 
    dplyr::filter(fcode == 46000 | fcode == 46003 | fcode == 46006 | fcode == 46007) %>% # feature codes for "stream/river"; see https://files.hawaii.gov/dbedt/op/gis/data/NHD%20Complete%20FCode%20Attribute%20Value%20List.pdf
    st_crop(box_name)
}
# Apply convert_nhd() to all oregon nhd filepaths at once
# new_nhd_or <- lapply(nhd_or_ls, convert_nhd, box_name = box_or)
# Combine resulting list of sf objects into a single dataframe/sf object
# nhd_or <- do.call(rbind, new_nhd_or)
# Save for later (above process take a long time)
#st_write(nhd_or, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/OR/nhd_or_coast.shp")

# Repeat this process for Washington state
# nhd_wa_ls <- intersect(list.files(path='/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/WA/Shape/',pattern='.shp',full.names = TRUE), list.files(path='/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/WA/Shape/', pattern="NHDFlowline", full.names = TRUE))
# box_wa <- c(xmin = -124.8712, ymin = 45.15649, xmax = -122, ymax = 50.62271)
# new_nhd_wa <- lapply(nhd_wa_ls, convert_nhd, box_name = box_wa)
# nhd_wa <- do.call(rbind, new_nhd_wa)

# And CA is already in sf object, which is:
# nhd_ca <- st_read("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/Rivers/CA/California_Streams.shp") %>% st_transform(crs = 4326)
# st_write(nhd_ca, "/Users/rachelcarlson/Desktop/foo1.shp")

# CALIFORNIA: Find nearest stream to each point after filtering points to only nearshore points (<= 10 m from shore)
full_daily_ca_coast <- full_daily_ca %>% dplyr::filter(distance_offshore <= 10)
nearstr_ca <- st_nearest_feature(full_daily_ca_coast, nhd_ca) # Find index value of nearest stream to each OA point
dist = st_distance(full_daily_ca_coast, nhd_ca[nearstr_ca,], by_element=TRUE) # Find distance between each OA point and nearest stream feature
# Combine OA dataset with nearest stream attributes
OA_ca <- cbind(full_daily_ca_coast, (nhd_ca[nearstr_ca,] %>% as.data.frame() %>% dplyr::select(c(DFGWATERID, Mouth_Lat, Mouth_Long, Down_Name, Down_ID, NHD_Perman, Hydrologic))))
# and distance to that stream
OA_ca$distance_stream <- dist

# The previous analysis just calculated distance to any point on linestring, but this might not be where river lets out.
# Try a separate metric for distance from nearest river mouth.
nhd_ca_mouth <- nhd_ca %>% as.data.frame() %>% dplyr::select(-geometry) %>% st_as_sf(coords = c("Mouth_Long", "Mouth_Lat"), crs = crs(full_daily_ca_coast))
nearmth_ca <- st_nearest_feature(full_daily_ca_coast, nhd_ca_mouth)
dist2 = st_distance(full_daily_ca_coast, nhd_ca_mouth[nearmth_ca,], by_element=TRUE)
OA_ca$distance_mouth <- dist2

OA_ca <- OA_ca %>% st_as_sf()

# ALL STATES: Combine all states' stream/river info into one dataframe. Stream/river attributes will be significantly pared down because CA, OR, WA datasets have inconsistent attributes. Use the CA-specific dataset above it need to access more granular CA stream/river information.
nhd_ca_trim <- nhd_ca %>% dplyr::select(c(NHD_Perman, SHAPESTLen, geometry))
nhd_ca_trim$fcode <- 46006
nhd_or_trim <- nhd_or %>% dplyr::select(c(permanent_, lengthkm, fcode, geometry))
nhd_wa_trim <- nhd_wa %>% dplyr::select(c(permanent_, lengthkm, fcode, geometry))
colnames(nhd_or_trim) <- c("NHD_Perman", "lengthkm", "fcode", "geometry")
colnames(nhd_wa_trim) <- c("NHD_Perman", "lengthkm", "fcode", "geometry")
colnames(nhd_ca_trim) <- c("NHD_Perman", "lengthkm", "geometry", "fcode")

ls_nhds <- list(nhd_ca_trim, nhd_wa_trim, nhd_or_trim)
nhd_all <- do.call(rbind, ls_nhds)

full_daily_coast <- full_daily %>% dplyr::filter(distance_offshore <=10) # Filter OA daily dataset to those points within 10 km of coast
nearstr <- st_nearest_feature(full_daily_coast, nhd_all)
dist = st_distance(full_daily_coast, nhd_all[nearstr,], by_element=TRUE) # Find distance between each OA point and nearest stream feature
OA_join <- cbind(full_daily_coast, (nhd_all[nearstr,] %>% as.data.frame()))
# and distance to that stream
OA_join$distance_stream <- dist


######## Calculate, for OA dataset, watershed-based alkalinity at HUC8, HUC10, and HUC12 levels
# Find nearest HUC8, 10, and 12 to each OA point
nearhuc8 <- st_nearest_feature(OA_join, huc8_ta)
# The avove is ONLY HUC8 that has nonzero alkalinity from inland data. For the flow attribute, you will also need nearest HUC8 regardless of whether/not there is alkalinity data.
nearhuc8.coast <- st_nearest_feature(OA_join, huc8)
nearhuc10 <- st_nearest_feature(OA_join, huc10_ta)
nearhuc12 <- st_nearest_feature(OA_join, huc12_ta)
# We don't need to know the distance to the watershed/HUC because we have distance to nearest stream
OA_join <- cbind(OA_join, (huc8_ta[nearhuc8,] %>% as.data.frame() %>% dplyr::select(-geometry)))
OA_join <- cbind(OA_join, (huc10_ta[nearhuc10,] %>% as.data.frame() %>% dplyr::select(-geometry)))
OA_join <- cbind(OA_join, (huc12_ta[nearhuc12,] %>% as.data.frame() %>% dplyr::select(-geometry)))
# Nearest coastal HUC8 (rather than nearest HUC8 with non-NA alkalinity)
index <- huc8[nearhuc8.coast,] %>% as.data.frame() %>% dplyr::select(-geometry)
OA_join$HUC8.coastal <- index$HUC8

######## Calculate flow rate per watershed on date of OA sampling
ca_flow <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/flow/All/CA_coast", header = FALSE, sep = "")
wa_flow <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/flow/All/WA_coast", header = FALSE, sep = "")
or_flow <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/flow/All/OR_coast", header = FALSE, sep = "")

# Remove extra text (before data) from text file
ca_flow <- ca_flow[-(1:288),]
wa_flow <- wa_flow[-c(1:44),]
or_flow <- or_flow[-(1:59),]

colnames(ca_flow) <- c("agency", "site.code", "date", "time", "timezone", "flow", "Qcode")
colnames(wa_flow) <- c("agency", "site.code", "date", "time", "timezone", "flow", "Qcode")
colnames(or_flow) <- c("agency", "site.code", "date", "time", "timezone", "flow", "Qcode")

ca_flow <- ca_flow[,-(8:15)]
wa_flow <- wa_flow[,-(8:15)]
or_flow <- or_flow[,-(8:15)]

ca_flow$date <- as.Date(ca_flow$date, format = "%Y-%m-%d")
wa_flow$date <- as.Date(wa_flow$date, format = "%Y-%m-%d")
or_flow$date <- as.Date(or_flow$date, format = "%Y-%m-%d")

ca_flow$flow <- as.numeric(ca_flow$flow)
wa_flow <- wa_flow %>% filter(!is.na(flow))
wa_flow$flow <- as.numeric(wa_flow$flow)
or_flow$flow <- as.numeric(or_flow$flow)

# Summarize watershed flow by date
ca_flow_daily <- ca_flow %>% group_by(site.code, date) %>% summarize(mnflow = mean(flow, na.rm = TRUE))
wa_flow_daily <- wa_flow %>% group_by(site.code, date) %>% summarize(mnflow = mean(flow, na.rm = TRUE))
or_flow_daily <- or_flow %>% group_by(site.code, date) %>% summarize(mnflow = mean(flow, na.rm = TRUE))

ca_flow_daily$month <- month(ca_flow_daily$date)
wa_flow_daily$month <- month(wa_flow_daily$date)
or_flow_daily$month <- month(or_flow_daily$date)

ca_flow_month.av <- ca_flow_daily %>% group_by(site.code, month) %>% summarize(mnflow = mean(mnflow, na.rm = TRUE),
                                                                               mxflow = max(mnflow, na.rm = TRUE))

# Bind flow to main dataset
ls_flow <- list(ca_flow_daily, wa_flow_daily, or_flow_daily)
flow_all <- do.call(rbind, ls_flow)
# Associate each site code with a HUC8
ca.sitecode <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/flow/All/CA_wqsite_by_huc.csv")
wa.sitecode <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/flow/All/WA_wqsite_by_huc.csv")
or.sitecode <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/Data/Hydrology/flow/All/OR_wqsite_by_huc.csv")

ls_sitecode <- list(ca.sitecode, wa.sitecode, or.sitecode)
sitecode_all <- do.call(rbind, ls_sitecode)
colnames(sitecode_all)[1] <- "site.code"

flow_all2 <- flow_all
flow_all2$site.code <- as.numeric(flow_all2$site.code)
flow_all3 <- flow_all2 %>% 
  right_join(sitecode_all, by = "site.code") %>% 
  group_by(HUC8, date) %>% 
  na.omit() %>% 
  summarize(mnflow = mean(mnflow, na.rm = TRUE),
            mxflow = max(mxflow, na.rm = TRUE))



foo <- OA_join %>% inner_join(flow_all, by=c('HUC8','date'))

# CALIFORNIA ONLY HUC DATA (IF NEED MORE GRANDULAR DATA)
# Find nearest HUC8, 10, and 12 to each OA point
# nearhuc8_ca <- st_nearest_feature(OA_ca, huc8_ta)
# nearhuc10_ca <- st_nearest_feature(OA_ca, huc10_ta)
# nearhuc12_ca <- st_nearest_feature(OA_ca, huc12_ta)
# # We don't need to know the distance to the watershed/HUC because we have distance to nearest stream
# OA_ca <- cbind(OA_ca, (huc8_ta[nearhuc8_ca,] %>% as.data.frame() %>% dplyr::select(-geometry)))
# OA_ca <- cbind(OA_ca, (huc10_ta[nearhuc10_ca,] %>% as.data.frame() %>% dplyr::select(-geometry)))
# OA_ca <- cbind(OA_ca, (huc12_ta[nearhuc12_ca,] %>% as.data.frame() %>% dplyr::select(-geometry)))
# 
# colnames(OA_ca[,c("mean_ta_umolkg",)])
# Filter out values where different parameter codes show radically different alkalinities (shown by too low or high of quotient between paramter values)
# new_usgs$screen1 <- new_usgs$p00418/new_usgs$p29801
# new_usgs$screen2 <- new_usgs$p29802/new_usgs$p29801
# new_usgs$screen3 <- new_usgs$p29802/new_usgs$p00418
# 
# new_usgs_ca <- new_usgs %>% filter((screen1 >= 0.75 & screen1 <= 1.25) | (screen2 >= 0.75 & screen2 <= 1.25) | (screen3 >= 0.75 & screen3 <= 1.25) | is.na(screen1))
# plot(new_usgs_ca$p00418, new_usgs_ca$p29802, xlim = c(0,2500), ylim = c(0,2500))
# 
# new_usgs_ca$alkalinity <- ifelse((new_usgs_ca$screen1 >= 0.75 & new_usgs_ca$screen1 <= 1.25), (new_usgs_ca$p00418 + new_usgs_ca$p29801)/2,
#                                  ifelse((new_usgs_ca$screen2 >= 0.75 & new_usgs_ca$screen2 <= 1.25), (new_usgs_ca$p29802 + new_usgs_ca$p29801)/2,
#                                         ifelse((new_usgs_ca$screen3 >= 0.75 & new_usgs_ca$screen3 <= 1.25), (new_usgs_ca$p29802 + new_usgs_ca$p00418)/2, NA)))
# 
# new_usgs_ca$alkalinity <- ifelse(is.na(new_usgs_ca$alkalinity), new_usgs_ca$p00418, new_usgs_ca$alkalinity)
# new_usgs_ca$alkalinity <- ifelse(is.na(new_usgs_ca$alkalinity), new_usgs_ca$p29801, new_usgs_ca$alkalinity)
# new_usgs_ca$alkalinity <- ifelse(is.na(new_usgs_ca$alkalinity), new_usgs_ca$p29802, new_usgs_ca$alkalinity)

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







# Import sites tagged for river proximity
# fresh <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/planning/sites.csv")
# fresh$Lat <- as.numeric(sub("0+$", "", as.character(fresh$Lat))) %>% round(digits = 4)
# fresh$Lon <- as.numeric(sub("0+$", "", as.character(fresh$Lon))) %>% round(digits = 4)
# fresh <- fresh %>% dplyr::select(c("Site", "Dataset_ID","Distance","Reach_ID","Flow_range_midFeb","Mean_flow","Lat","Lon"))
# full <- full %>% rename("Lat" = "latitude", "Lon" = "longitude", "Dataset_ID" = "dataset_id") # There were 70ish mismatched dataset IDs in first try at joining these datasets. Ensure datest ID in tagged freshwater dataset is correct.

# Join reach, site name, flow, and distance information to relevant sites
# full_fresh <- full %>% full_join(fresh, by = c("Lat" = "Lat", "Lon" = "Lon", "Dataset_ID" = "Dataset_ID"), multiple = "all")
# full_fresh$Distance <- as.numeric(full_fresh$Distance)
# Join additional Columbia sites from WCOA cruises. Assign "Columbia River" as site name so that these cruise data are included.
# full_fresh$Site <- ifelse((full_fresh$Lat >= 46) & (full_fresh$Lat <= 46.3) & (full_fresh$Lon >= -124.22) & (full_fresh$depth_m <= 10) & (!is.na(full_fresh$pH_total)), "Columbia River", full_fresh$Site)
# full_fresh$distance_offshore <- ifelse((full_fresh$Lat >= 46) & (full_fresh$Lat <= 46.3) & (full_fresh$Lon >= -124.22) & (full_fresh$Dataset_ID == 21) & (!is.na(full_fresh$pH_total)), "", full_fresh$distance_offshore)
### List to search inland
# Goleta is ideal
# Cambria/Hearst State Park (Santa Rosa Creek and Leffingwell Creek; Crystal Cove (2); San Clamente

full_fresh$date <- as.Date(full_fresh$time_utc)

cond_fresh <- full_fresh %>% filter(!is.na(ta_umolkg)) %>% filter(!is.na(Site)) %>% 
  group_by(date, Site, Dataset_ID, Lat, Lon) %>% summarize(pH_daily = mean(pH_total, na.rm = TRUE),
                                                           temp_daily = mean(t_C, na.rm = TRUE),
                                                           alk_daily = mean(ta_umolkg, na.rm = TRUE),
                                                           sal_daily = mean(sal_pss, na.rm = TRUE),
                                                           depth_mean = mean(depth_m, na.rm = TRUE),
                                                           distance = mean(Distance, na.rm = TRUE))

write.csv(cond_fresh, "/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/OA_all_alkalinity.csv")


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


Columbia1 <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/planning/66_Final_OOI_SLHNov2022.csv")
OOI <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/planning/66_Final_OOI_SLHNov2022.csv") %>% dplyr::select(-X)# OOI data (dataset 66)
Trinidad <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/planning/67_TrinidadHeadLine/67_TrinidadHead_Nov2022SLH.csv") %>% dplyr::select(-X)
Newport <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/planning/68_NewportHydrographicLine/68_NewportHydrographicLine_Final.csv") %>% dplyr::select(-c(X,X.1))
Columbia <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/planning/70_CMOPSaturn2/70_CMOPSaturn2_final.csv") %>% dplyr::select(-X)
Scripps <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/OA_dataset/planning/71_ScrippsLaJolla/71_ScrippsLaJolla_final.csv") %>% dplyr::select(-X)
variable.names(full)
variable.names(Scripps)



# 5, 23, 21

# Edit May 2: I don't think the below is actually true:
# Dataset 27 pH values are filtered out of the 'clean' dataset because they use YSI. However, this dataset will be important to us.
# Filter out "clean" version of dataset 27 and then load the full (unfiltered) version back in.
# full <- full %>% filter(dataset_id != 27)
# CoastOA <- read.csv("/Users/rachelcarlson/Documents/Research/Postdoc-2022-present/Freshwater_time/planning/27_CoastalOA/27_final_coastal_OA.csv")

# Bind datasets
# OOI$distance_offshore <- NA # Create a blank column for distance_offshore as this is missing in individual datasets
# Trinidad$distance_offshore <- NA
# Newport$distance_offshore <- NA
# Columbia$distance_offshore <- NA
# Scripps$distance_offshore <- NA
# CoastOA$distance_offshore <- NA

# full2 <- do.call("rbind", list(full, OOI, Trinidad, Newport, Columbia, Scripps, CoastOA)) # Bind datasets

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

