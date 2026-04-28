######### Exploratory Data Analysis (EDA)
######### Freshwater Influence on Ocean Acidification
######### 27 November 2026
######### Rachel R. Carlson

# ==============================================================================
# SECTION 1: DATA PREPARATION AND EXPLORATION
# ==============================================================================

# -----------------------------------------------------------------------------
# 1.1 Load Required Libraries
# -----------------------------------------------------------------------------
library(tidyverse)
library(ncdf4)
library(raster)
library(sf)
library(purrr)
library(lubridate)
library(lme4)
library(foreign) # allows you to read DBF
library(readr)
library(corrplot)
library(lattice)
library(mgcv)
library(scales)
library(ggplot2)
library(spdep)    # For spatial autocorrelation tests
library(gstat)    # For variogram analysis
library(lmtest)   # For heteroscedasticity tests

# -----------------------------------------------------------------------------
# 1.2 Load Dataset
# -----------------------------------------------------------------------------
# Read the clustered OA dataset from Data_wrangling.R output
oa <- read.csv("/Users/rachelcarlson/Documents/Berkeley/Research/Postdoc-2022-present/OA_dataset/Data/OA/OA_master_clustered_26April2026.csv") %>% dplyr::select(-X)

# Convert distance from meters to kilometers (fixing unit mistake in previous code)
oa$distanceRiv_km <- oa$distanceRiv / 1000

# -----------------------------------------------------------------------------
# 1.3 Create State Variable for Regional Analysis
# -----------------------------------------------------------------------------
# Categorize observations by US West Coast state based on latitude
oa <- oa %>%
  mutate(state = case_when(
    lat < 42 ~ "CA",           # California: lat < 42
    lat >= 42 & lat < 46 ~ "OR", # Oregon: 42 <= lat < 46
    lat >= 46 ~ "WA",           # Washington: lat >= 46
    TRUE ~ NA_character_
  ))

# -----------------------------------------------------------------------------
# 1.4 Correlation Matrix: Assess Collinearity Among Covariates
# -----------------------------------------------------------------------------
# Select numeric independent variables for correlation analysis
# These represent potential predictors for OA metrics
ColMatrix <- oa %>% dplyr::select(c(seaSal, distanceRiv, composite_fw_expose, 
                                     soilpH_hucMax, depth_bin, upwellN, upwell, 
                                     month, year))

# Display structure of selected variables
str(ColMatrix)

# Compute Pearson correlation matrix (pairwise complete observations)
ColMatrix_cor <- cor(ColMatrix, use = "pairwise.complete.obs", method = "pearson")

# Visualize correlation matrix with colored tiles and numeric coefficients
# Note: Covariates with correlation >= 0.7 should have only one retained 
# to avoid multicollinearity (e.g., max and mean HUC8 total alkalinity)
corrplot(ColMatrix_cor, method = "color", type = "upper", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45, number.cex = 0.7)

# -----------------------------------------------------------------------------
# 1.5 Exploratory Plots: Freshwater Influence on Ocean Chemistry
# -----------------------------------------------------------------------------

# Plot 1: Total Alkalinity vs. Ocean Salinity
# Freshwater has lower salinity, so this examines how freshwater influence 
# affects alkalinity in coastal ocean waters
sum(is.na(oa$seaSal)) # Check for missing salinity values (37 missing)

ggplot(oa, aes(x = seaSal, y = seaTA)) +
  geom_point(color = "#474540", size = 1.5) +
  geom_smooth(method = "lm") +
  theme_light() + 
  xlab("Salinity (ppt)") + ylab("Total Alkalinity (umol kg-1)")

# Plot 2: TA vs. Salinity colored by freshwater exposure index
ggplot(oa, aes(x = seaSal, y = seaTA, color = log(composite_fw_expose))) +
  geom_point(size = 1.5) +
  scale_color_viridis_c(name = "FW exposure index", option = "viridis", direction = -1) +
  geom_smooth(method = "lm") +
  theme_light() + 
  xlab("Salinity (ppt)") + ylab("Total Alkalinity (umol kg-1)")

# Plot 3: TA vs. Salinity colored by soil pH from watershed. Soil pH influences alkalinity through watershed drainage
# Shows some signs of negative residuals (lower than expected TA) in areas with low soil pH.
ggplot(oa, aes(x = seaSal, y = seaTA, color = soilpH_hucMean)) +
  geom_point(size = 1.5) +
  scale_color_viridis_c(name = "Soil pH in Water", option = "viridis", direction = -1) +
  geom_smooth(method = "lm") +
  theme_light() + 
  xlab("Salinity (ppt)") + ylab("Total Alkalinity (umol kg-1)")

# Plot 4: Same as plot 3, but filter to mixed layer only (depth_bin = 0, i.e., 0-15 m)
ggplot(oa %>% filter(depth_bin == 0), aes(x = seaSal, y = seaTA, color = soilpH_hucMean)) +
  geom_point(size = 1.5) +
  scale_color_viridis_c(name = "Soil pH in Water", option = "viridis", direction = -1) +
  geom_smooth(method = "lm") +
  theme_light() + 
  xlab("Salinity (ppt)") + ylab("Total Alkalinity (umol kg-1)")


# -----------------------------------------------------------------------------
# # 1.7 Additional Exploratory Plots
# # -----------------------------------------------------------------------------
# 
# # Plot: Site-level TA vs. Salinity colored by distance from river
# ggplot(oa %>% filter(distanceRiv <= 500), aes(x = seaSal, y = seaTA, color = distanceRiv)) +
#   geom_point(size = 1.5) +
#   geom_smooth(method = "lm", color = "black") +
#   scale_color_gradient(low = "#2c7fb8", high = "#f03b20", name = "Distance from River") +
#   theme_light() + 
#   xlab("Salinity (ppt)") + ylab("Total Alkalinity (umol kg-1)")
# 
# # Plot: Freshwater exposure index vs. Ocean Salinity
# ggplot(oa %>% filter(distanceRiv <= 500), aes(x = composite_fw_expose_month, y = seaSal, color = seaTA)) +
#   geom_point(size = 1.5) +
#   scale_color_gradient(low = "#2c7fb8", high = "#f03b20", name = "Alkalinity") +
#   theme_light() + 
#   xlab("Freshwater exposure index") + ylab("Salinity (ppt)")
# 
# # Plot: Distance from shore vs. Salinity
# ggplot(shallow_siteMeans, aes(x = distanceRiv, y = seaSal, color = distanceRiv_km)) +
#   geom_point(size = 1.5) +
#   scale_color_gradient(low = "#2c7fb8", high = "#f03b20", name = "Distance from River") +
#   theme_light() + 
#   xlab("Distance from Shore (m)") + ylab("Salinity (ppt)")



# ==============================================================================
# SECTION 2: STATISTICAL MODELING
# ==============================================================================

# -----------------------------------------------------------------------------
# 2.1 Data Rescaling for Regression
# -----------------------------------------------------------------------------
# Filter to mixed layer only (depth_bin == 0) and rescale covariates
# Rescaling improves model convergence and interpretability

oas <- oa %>% filter(depth_bin == 0)
oas$depth_bin <- rescale(log(oas$depth_bin + 0.0000001))
oas$composite_fw_expose <- rescale(log(oas$composite_fw_expose + 0.0000001))
oas$soilpH_hucMean <- rescale(oas$soilpH_hucMean)
oas$soilpH_hucMax <- rescale(oas$soilpH_hucMax)
oas$upwell <- rescale(oas$upwell)
oas$upwellN <- rescale(oas$upwellN)
oas$lat <- rescale(oas$lat)
oas$seaSal <- rescale(oas$seaSal)
oas$distanceRiv <- rescale(oas$distanceRiv)
oas$distanceSho <- rescale(sqrt(oas$distanceSho))

# -----------------------------------------------------------------------------
# 1.8 Variable Distribution Checks
# -----------------------------------------------------------------------------
# Check distributions of key variables in the filtered dataset (mixed layer only)

# Boxplot: TA by dataset ID (to identify dataset-specific patterns)
# Dataset 27 and 38 have by far the most variable TA
boxplot(seaTA~dataset_id, data=oas, main="TA by dataset ID",
     xlab="Dataset ID", ylab="TA") 

# Scatterplots: TA vs. distance variables
# As expected, TA is most variable close to shore
plot(seaTA~distanceRiv, data=oas, main="TA by distance from river",
     xlab="Distance from River", ylab="TA")

# Also as expected, TA is most variable in the shallows (freshwater floats due to low density)
plot(seaTA~depth_m, data=oas, main="TA by depth",
     xlab="Depth (m)", ylab="TA")

# Boxplot: TA by month (seasonal patterns)
oas$fMonth <- factor(oas$month)
bwplot(seaTA~fMonth | dataset_id, data=oas, xlab = "Month")

# Boxplot: TA by year (temporal trends)
oas$fYear <- factor(oas$year)
bwplot(seaTA~fYear | dataset_id, data=oas, xlab = "Year")

# -----------------------------------------------------------------------------
# 2.2 Linear Model (LM)
# -----------------------------------------------------------------------------
# Basic linear model examining relationships between ocean TA and predictors
# Includes interaction between salinity and soil pH (watershed influence)

foo <- lm(seaTA ~ seaSal*soilpH_hucMax + upwell + composite_fw_expose + state + month + year, data = oas)
plot(foo)

# My hunch is that TA and Salinity are positively related (this is a well known relationship) but that there is lower than expected SeaTA near watersheds with low soil pH. The last part would be our main finding.
# However, the diagnostics for all of these models are terrible.

# -----------------------------------------------------------------------------
# 2.3 Generalized Linear Model (GLM) with Gamma Distribution
# -----------------------------------------------------------------------------
# Gamma regression is appropriate for positive continuous response variables
# with right-skewed distributions (like Total Alkalinity)

# Fit gamma GLM with log link (common choice for gamma models)
glm.gamma <- glm(seaTA ~ seaSal*soilpH_hucMax + upwell + composite_fw_expose + state + depth_bin,
                 data = oas, 
                 family = Gamma(link = "log"))

# Summary of the model
summary(glm.gamma)

# Alternative: gamma with inverse link (for comparison)
glm.gamma.inv <- glm(seaTA ~ seaSal*soilpH_hucMax + upwell + composite_fw_expose + state + depth_bin,
                     data = oas, 
                     family = Gamma(link = "inverse"))

summary(glm.gamma.inv)

# -----------------------------------------------------------------------------
# 2.4 Diagnostic Tests for Gamma GLM
# -----------------------------------------------------------------------------

# Test 1: Check for overdispersion
# Gamma models can exhibit overdispersion; check the dispersion parameter
dispersion <- summary(glm.gamma)$dispersion
cat("Dispersion parameter:", dispersion, "\n")
# Dispersion is not >> 1, no need for quasi-poisson or negative binomial

# Test 2: Residual diagnostics
# Pearson residuals
pearson_resid <- residuals(glm.gamma, type = "pearson")
# Deviance residuals
deviance_resid <- residuals(glm.gamma, type = "deviance")

# Plot residuals vs fitted values
plot(fitted(glm.gamma), pearson_resid,
     xlab = "Fitted values", ylab = "Pearson residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2) # Residuals are tapering at 2000 - 2200, as in lm. 
# There's probably lots of spatial autocorrelation and nonlinear relationships (gam needed?)
# Another issue: some sites are revisited several times and others are not.

# Q-Q plot for residuals (it looks terrible)
qqnorm(deviance_resid, main = "Q-Q plot of deviance residuals")
qqline(deviance_resid)

# Test 3: Test for spatial autocorrelation in residuals
# Convert to sf object for spatial analysis
oas_sf <- st_as_sf(oas, coords = c("lon", "lat"), crs = 4326)

# Moran's I using spdep package
# Get rows used in the model
used_rows <- as.numeric(rownames(model.frame(glm.gamma)))
# Subset sf object to model rows
oas_sf_mod <- oas_sf[used_rows, ]
# Residuals from same model
deviance_resid <- residuals(glm.gamma, type = "deviance")
# Coordinates for only model rows
coords_nbd <- st_coordinates(oas_sf_mod)
# Find neighbors within 10km (10000m)
nb <- dnearneigh(coords_nbd, d1 = 0, d2 = 10000)
# Create spatial weights matrix
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
# Transform to projected CRS (meters; good for CONUS)
oas_sf <- st_transform(oas_sf, 5070)
# Moran's I test
moran.test(deviance_resid, lw, zero.policy = TRUE)

# Test 4: Variogram at different lags (alternative spatial test)
variogram_resid <- variogram(deviance_resid ~ 1, locations = ~longitude + latitude, data = oas)
plot(variogram_resid, main = "Variogram of residuals")

# Test 5: Test for heteroscedasticity (Breusch-Pagan test)
bptest(glm.gamma)

# Test 6: Component + residual plots (test for linearity)
crPlots(glm.gamma)

# Test 7: Influence measures (Cook's distance)
plot(cooks.distance(glm.gamma), type = "h",
     main = "Cook's Distance", ylab = "Cook's D")
abline(h = 4/nrow(oas), lty = 2)  # Common threshold

# Test 8: Compare models using AIC/BIC
AIC(glm.gamma, glm.gamma.inv)
BIC(glm.gamma, glm.gamma.inv)

# -----------------------------------------------------------------------------
# 2.5 Generalized Additive Models (GAM)
# -----------------------------------------------------------------------------
# GAMs allow for non-linear relationships between predictors and response
# Using the same variable names as the LM and GLM models above

# GAM model 1: Total alkalinity as response (using same variables as LM/GLM)
gam.m1 <- gam(seaTA ~ s(distanceRiv) + s(composite_fw_expose) + s(soilpH_hucMax) + 
              s(upwell) + s(month) + s(year) + state,
              data = oas, family = gaussian, method = "REML")

# GAM model 2: Alternative specification with reduced smooth terms
gam.m2 <- gam(seaTA ~ s(seaSal) + s(soilpH_hucMax) + s(upwell) + s(composite_fw_expose) + 
              state + month + year,
              data = oas, family = gaussian, method = "REML")

# GAM model 3: Interaction between salinity and soil pH (by state)
gam.m3 <- gam(seaTA ~ te(seaSal, soilpH_hucMax) + s(composite_fw_expose) + 
              upwell + state + s(month, bs = "cc") + s(year),
              data = oas, family = gaussian, method = "REML")

# -----------------------------------------------------------------------------
# 2.6 Diagnostic tests for Generalized Additive Models (GAM)
# -----------------------------------------------------------------------------












