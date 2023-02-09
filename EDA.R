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
id_list <- c(30, 31, 32, 33, 39, 40, 50, 56, 64, 66, 67, 70)
trim <- full %>% filter(dataset_id %in% id_list)
trim <- trim %>% group_by(dataset_id) %>% mutate(min_time = min(time_utc, na.rm = TRUE),
                                            max_time = max(time_utc, na.rm = TRUE))

View(trim %>% group_by(dataset_id) %>% summarize(min_time = min(time_utc, na.rm = TRUE),
                                         max_time = max(time_utc, na.rm = TRUE)))

# map the sites
coords1 <- trim %>% distinct(latitude, longitude, dataset_id, min_time, max_time)
mcoords <- st_as_sf(coords1, coords = c("longitude", "latitude"), crs = 4326)
plot(mcoords$geometry)

# Calculate distance from shore. Sf isn't working here, so exported to ArcGIS where I ran the "near" tool on all sites.



