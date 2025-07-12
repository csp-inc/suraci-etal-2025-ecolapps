## ---------------------------
##
## Script name: 02-prep-ssf-dataset.R
##
## Purpose: Perpare analysis-ready used-
## available datasets for SSF models
##
## Author: Justin Suraci
##
## Email contact: justin@csp-inc.org
##
## ---------------------------

library(tidyverse)
library(sf)
library(terra)
library(amt)
library(lubridate)

# ---------- Collect Datasets -------------------
# Get AOI
puma_crs <- st_crs(26910)
aoi <- st_read("data/puma-refined-aoi.geojson") %>% st_transform(puma_crs)

# Get puma location data
load(file = 'data/puma-locs-full-age-class-20230321.rda', verbose = TRUE)
#subset to AOI
puma_sp <- st_as_sf(puma_ac, coords = c("UTM_X", "UTM_Y"), crs = puma_crs)
puma_in <- st_within(puma_sp, aoi, sparse = F) %>% as.vector()
puma_ac <- puma_ac[puma_in,]
# identify individuals with relocation time <= 2 hrs
reloc <- puma_ac %>% group_by(name) %>% 
  summarize(hr_median = (median(dt_hr))) %>%
  filter(hr_median <=2.5) %>% as.data.frame()

# Process puma location dataset
puma_locs <- puma_ac %>%
  # Remove dependent young
  filter(age_class != "dep_young",
         name %in% reloc$name) %>% 
  # Create age + year + id column to break tracks up by age class and year (for yearly covs)
  mutate(yr = year(date_time),
         age_ID = paste(name, age_class, yr, sep = "_")) 

# Get and organize raster covs
dcovs <- rast("data/ssf-cov-rasters/puma-ssf-dist-covs-270m.tif") # Import 270m covariates
lcovs <- rast("data/ssf-cov-rasters/puma-ssf-focal-covs-30m.tif") # Import 30m covariates
fcovs <- rast("data/ssf-cov-rasters/puma-ssf-forestry-covs-30m.tif")
# Subset covs for extraction
covsPt_fine <- lcovs[[c("elevation_150m", "elevation_1km", "shrub_pcov_150m", "shrub_pcov_1km",
                     "ag_pcov_150m", "ag_pcov_1km", "developed_150m", "developed_1km",
                     "disturbed_edge_pcov_150m", "disturbed_edge_pcov_1km", "human_pop_den_150m",
                     "human_pop_den_1km", "nightlight_150m", "nightlight_1km")]]
covsPt_fine <- c(covsPt_fine, fcovs)
covsPt_course <- dcovs[[c("riparian_dist", "ag_dist", "road_dist", "developed_dist")]]
covsPt_year <- lcovs[[str_detect(names(lcovs), "20..edge_pcov|20.._trCover_pcov")]]
covsPath <- lcovs[[c("slope_150m", "vrm", "tpi", "impervious_pcov_150m", 
                     "developed_150m", "ag_pcov_150m")]]
covsPath <- c(covsPath, fcovs["forestry_pcov_150m"])


# ---------- Extract Covs -------------------
path_extractor = function(steps){
  coNames <- paste(names(steps), "path", sep = "_")
  coMeans <- apply(steps, 2, mean, na.rm = T)
  names(coMeans) <- coNames
  return(coMeans)
}

ids <- unique(puma_locs$age_ID)
ssf_ua <- list()
for(i in 1:length(ids)){
  
  print(paste("running animal-year #", i, "of", length(ids)))
  
  # Subset to focal animal and add "previous step length" to remove short (likely stationary) steps
  temp <- puma_locs %>% filter(age_ID == ids[i])
  temp$pstep <- c(NA, temp$step)[1:nrow(temp)]
  temp <- temp[which(is.na(temp$pstep) | temp$pstep > 20),]
  
  # Get appropriate yearly covs and remove year from name
  yr <- temp$yr[1] %>% as.character()
  covsYr <- covsPt_year[[str_detect(names(covsPt_year), yr)]] 
  cnames <- names(covsYr) %>% str_split(yr) %>% map_chr(2) %>% str_replace("_tr", "tr")
  names(covsYr) <- cnames
  
  # Create used and available steps and extract raster covs at end points of each
  ua_temp <- temp %>% make_track(UTM_X, UTM_Y, date_time, name = name, sex = sex, yr = yr,
                                  age_class = age_class, moon_frac = moon_frac, moon_alt = moon_alt,
                                  sun_alt = sun_alt, day_night = day_night, locID = locID, crs = 26910) %>% 
    track_resample(rate = hours(2), tolerance = minutes(30)) %>% 
    steps_by_burst(keep_cols = 'end')
  # Skip puma-years with too few steps
  if(dim(ua_temp)[1]<30) next
  ua_temp <- ua_temp %>% 
    random_steps() %>% 
    extract_covariates(covsPt_fine, where = 'end') %>%
    extract_covariates(covsPt_course, where = 'end') %>%
    extract_covariates(covsYr, where = 'end') %>%
    mutate(cos_ta = cos(ta_),
           log_sl = log(sl_))
  
  # Get along path covariates and bind to above dataset
  path_extract <- ua_temp %>% extract_covariates_along(covsPath) %>% map(path_extractor) %>% do.call(rbind,.)
  ua_temp <- ua_temp %>% cbind(path_extract)
  
  # add to list
  ssf_ua[[i]] <- ua_temp
}
puma_ssf_dat <- do.call(rbind, ssf_ua)

# Add fully unique step id
puma_ssf_dat <- puma_ssf_dat %>% 
  mutate(unique_step = paste(name, age_class, yr, step_id_, sep = "_"))

# Save it
save(puma_ssf_dat, file = "data/puma-ssf-dataset-final.rda")
write_csv(puma_ssf_dat, file = "data/puma-ssf-dataset-final.csv")


