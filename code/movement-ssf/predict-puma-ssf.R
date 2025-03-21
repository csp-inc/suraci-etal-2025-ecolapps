library(tidyverse)
library(sf)
library(terra)
library(amt)
library(lubridate)
library(glmmTMB)

#------------- Compile data/models ------------------

# Get AOI
puma_crs <- st_crs(26910)
aoi = st_read("data/puma-refined-aoi.geojson") %>% st_transform(puma_crs)

# Load top models for each age/sex class
load("output/ad-fem-ssf-top-mod-v2.rda")
load("output/ad-male-ssf-top-mod-v2.rda")
load("output/disp-fem-ssf-top-mod-v2.rda")
load("output/disp-male-ssf-top-mod-v2.rda")
load("output/all-cats-ssf-top-mod.rda")

# Analysis Dataset
load(file = "data/puma-ssf-dataset-20230701.rda", verbose = TRUE)
ssf_dat <- puma_ssf_dat %>% 
  dplyr::mutate(case = as.numeric(case_))

# Covariate rasters
# Read in
lcovs <- rast("data/ssf-cov-rasters/puma-ssf-focal-covs-30m.tif") # Import 30m covariates
dcovs <- rast("data/ssf-cov-rasters/puma-ssf-dist-covs-270m.tif") # Import 270m covariates
dc_resamp <- dcovs[[c("riparian_dist", "road_dist")]] %>% disagg(fact = 9) %>% crop(lcovs$slope_150m) # Resample course covs to 30m
# Compile
covs <- lcovs[[c("shrub_pcov_150m", "ag_pcov_150m", "developed_150m",
                 "slope_150m","2022_trCover_pcov_1km","2022edge_pcov_150m")]]
covs <- c(covs, dc_resamp)
names(covs) <- c("shrub_pcov_150m", "ag_pcov_150m", "developed_150m",
                 "slope_150m_path","trCover_pcov_1km","edge_pcov_150m",
                 "riparian_dist","road_dist")

#---------------- Functions ---------------------
# Prep age/sex class specific datasets
datPrep <- function(dat){
  dout <- dat %>% dplyr::select(case, name, sex, day_night, locID, unique_step, cos_ta, log_sl,
                                moon_frac, shrub_pcov_150m, ag_pcov_150m, developed_150m, 
                                disturbed_edge_pcov_150m, nightlight_150m, forestry_pcov_1km, 
                                riparian_dist, ag_dist, road_dist, developed_dist, edge_pcov_150m, 
                                trCover_pcov_1km, slope_150m_path, tpi_path) %>% 
    mutate(moon = scale(moon_frac),
           night = ifelse(day_night=='Night',1,0),
           shrub = scale(shrub_pcov_150m),
           shrub_sq = shrub^2,
           ag = scale(ag_pcov_150m),
           ag_sq = ag^2,
           dev = scale(developed_150m),
           dev_sq = dev^2,
           disturbed = scale(disturbed_edge_pcov_150m),
           nightlight = scale(nightlight_150m),
           forestry = scale(forestry_pcov_1km),
           d_riparian = scale(riparian_dist),
           d_ag = scale(ag_dist),
           d_ag_sq = d_ag^2,
           d_road = scale(road_dist),
           d_road_sq = d_road^2,
           d_dev = scale(developed_dist),
           d_dev_sq = d_dev^2,
           edge = scale(edge_pcov_150m),
           tree = scale(trCover_pcov_1km),
           tree_sq = tree^2,
           slope = scale(slope_150m_path),
           slope_sq = slope^2,
           tpi = scale(tpi_path))
  return(dout)
}

# Function to extract covariate values at all cell locations in study area
gridCovs <- function(covs, st_area, dat){
  
  # Get cell numbers for raster cells overlapping study area
  pred_cell_nums <- terra::cells(covs, vect(st_area))[,2]
  # Extract covariate values from across the covs raster for use in pred surface
  pred_cov_temp <- terra::extract(covs, pred_cell_nums) 
  pred_xy <- terra::xyFromCell(covs, pred_cell_nums) %>% 
    as.data.frame() %>% 
    rename(utm_x = x, utm_y = y) 
  pred_covs <- cbind(pred_cov_temp, pred_xy) %>% 
    drop_na()
  
  # Rescale prediction covs based on range of values in rsf dataset
  covScale <- function(pc = pred_covs, in_dat = dat, cname){
    return((pc[,cname] - mean(in_dat[,cname], na.rm = TRUE))/sd(in_dat[,cname], na.rm = TRUE))
  }
  
  pred_covs_scale <- data.frame(slope = covScale(cname = "slope_150m_path"),
                                dev = covScale(cname = "developed_150m"),
                                ag = covScale(cname = "ag_pcov_150m"),
                                shrub = covScale(cname = "shrub_pcov_150m"),
                                edge = covScale(cname = "edge_pcov_150m"),
                                tree = covScale(cname = "trCover_pcov_1km"),
                                d_road = covScale(cname = "road_dist"),
                                d_riparian = covScale(cname = "riparian_dist"),
                                log_sl = mean(dat$log_sl, na.rm = TRUE),
                                cos_ta = mean(dat$cos_ta, na.rm = TRUE)) %>% 
    mutate(slope_sq = slope^2, 
           dev_sq = dev^2,
           ag_sq = ag^2,
           tree_sq = tree^2,
           shrub_sq = shrub^2,
           d_road_sq = d_road^2,
           "dev:slope" = dev * slope,
           utm_x = pred_covs$utm_x,
           utm_y = pred_covs$utm_y)
  
  # Return scaled prediction covs
  return(pred_covs_scale)
}

# Function to predict from model equation at all cell locations in study area
gridPredict <- function(gcovs, mod, scale = "exp", constrain = NULL){
  coefs <- fixef(mod)$cond # get model coefficients
  gcovs_mat <- gcovs %>% dplyr::select(names(coefs)) %>% data.matrix() # grab applicable columns from grid covs
  # Create linear prediction (and transform, as specified in function call)
  pred <- (gcovs_mat %*% coefs)
  if(!is.null(constrain)) pred <- pConstrain(pred, constrain)
  if(scale == "logit") pred <- pred
  if(scale == "exp") pred <- exp(pred)
  if(scale == "response") pred <- exp(pred)/(1+exp(pred))
  pred_out <- data.frame(prediction = pred, utm_x = gcovs$utm_x, utm_y = gcovs$utm_y)
  return(pred_out)
}

# Function to rasterize predictions
gridRast <- function(dat, filename=NULL, scale = 90){
  # Convert prediction dataset to sf and then terra SpatVector
  pred_sp <- dat %>% 
    st_as_sf(coords = c("utm_x", "utm_y"), crs = puma_crs) 
  
  # # Create template raster
  r_temp <- rast()
  terra::ext(r_temp) <- terra::ext(pred_sp)
  terra::res(r_temp) <- scale
  terra::crs(r_temp) <- puma_crs$wkt
  
  # Rasterize
  pred_sp <- vect(pred_sp)
  
  if(is.null(filename)) {
    pred_rast <- terra::rasterize(pred_sp, covs[[1]], field = 'prediction', fun = 'mean') %>% 
      crop(r_temp)
  } else {
    pred_rast <- terra::rasterize(pred_sp, covs[[1]], field = 'prediction', fun = 'mean', 
                                  filename = filename, overwrite = TRUE) %>% crop(r_temp)
  } 
  
  return(pred_rast)
}

# Set scale
pred_scale = 30

# Adult Female
ad_fem <- ssf_dat %>% filter(age_class == "adult", sex == "Female") %>% datPrep()
ad_fem_covs <- gridCovs(covs, aoi, ad_fem)
ad_fem_pred <- gridPredict(ad_fem_covs, ad_fem_top)
ad_fem_path <- paste0('output/ad-fem-ssf-prediction-', pred_scale, 'm.tif')
ad_fem_rast <- gridRast(ad_fem_pred, filename = ad_fem_path, scale = pred_scale)

# Adult Male
ad_male <- ssf_dat %>% filter(age_class == "adult", sex == "Male") %>% datPrep()
ad_male_covs <- gridCovs(covs, aoi, ad_male)
ad_male_pred <- gridPredict(ad_male_covs, ad_male_top)
ad_male_path <- paste0('output/ad-male-ssf-prediction-', pred_scale, 'm.tif')
ad_male_rast <- gridRast(ad_male_pred, filename = ad_male_path, scale = pred_scale)

# Disperser Female
disp_fem <- ssf_dat %>% filter(age_class == "disperser", sex == "Female") %>% datPrep()
disp_fem_covs <- gridCovs(covs, aoi, disp_fem)
disp_fem_pred <- gridPredict(disp_fem_covs, disp_fem_top)
disp_fem_path <- paste0('output/disp-fem-ssf-prediction-', pred_scale, 'm.tif')
disp_fem_rast <- gridRast(disp_fem_pred, filename = disp_fem_path, scale = pred_scale)

# Disperser Male
disp_male <- ssf_dat %>% filter(age_class == "disperser", sex == "Male") %>% datPrep()
disp_male_covs <- gridCovs(covs, aoi, disp_male)
disp_male_pred <- gridPredict(disp_male_covs, disp_male_top)
disp_male_path <- paste0('output/disp-male-ssf-prediction-', pred_scale, 'm.tif')
disp_male_rast <- gridRast(disp_male_pred, filename = disp_male_path, scale = pred_scale)

# All Cats
all_cats <- ssf_dat %>% datPrep()
all_cats_covs <- gridCovs(covs, aoi, all_cats)
all_cats_pred <- gridPredict(all_cats_covs, all_cats_top)
all_cats_path <- paste0('output/all-cats-ssf-prediction-', pred_scale, 'm.tif')
all_cats_rast <- gridRast(all_cats_pred, filename = all_cats_path, scale = pred_scale)
