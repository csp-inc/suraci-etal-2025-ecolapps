library(amt)
library(tidyverse)
library(terra)
library(blme)
library(sf)
library(glmmTMB)
library(car)
library(AICcmodavg)
library(rgdal)
library(raster)
library(ggpubr)


# ------------------------------------------------------------------------------
# DATA PREP
# ------------------------------------------------------------------------------
# --- AREA OF INTEREST ---
# Define AOI and CRS
su_crs <- st_crs(26910)
aoi <- st_read("final-covars-20230920/model-predictions/spatial-data/refined_AOI_03Feb2023_CLIPtoLAND_21Sept2023.shp")
aoi <- st_set_crs(aoi[1], "EPSG:4326") # defining original CRS
aoi <- st_transform(aoi, CRS("+init=epsg:26910 +proj=utm +ellps=WGS84 +datum=WGS84 +no_defs")) # transforming to match raster

# --- KILL SITES --- 
# Get age/sex kill site datasets (for calculating covariate means and SDs)
load("age-sex-class-RSFs/age-sex-pooled-RSF/age_sex_pooled_RSF_covariate_extract_all.rda") #load("age-sex-class-RSFs/age-sex-pooled-RSF/age_sex_pooled_RSF_covariate_extract.rda") <- commented out version was cleaned and scaled for modeling
# remove unused covariates
selected_df <- selected_df %>% 
  select(-c(aspect_100m, elevation_100m, disturbed_pcov_1km, grass_pcov_1km))

# --- COVARIATES ---
forestCovars <- raster::stack("final-covars/puma-kill-site-covs-forestvars-30m.tif")
forestedgeCovars <- raster::stack("final-covars/puma-kill-site-covs-forestedgevars-30m.tif")
nsResamp <- raster::stack("final-covars/nsResamp.grd")
distResamp <- raster::stack("final-covars/distResamp.grd")

# only select covariates that are retained in the analysis
nonSeasonalCovars <- nsResamp[[c(2,17,28,29)]]
forestry <- forestCovars[[c(6,25)]]
names(forestry[[1]]) <- "RAP_treecover" # currently selecting 2022
forestEdge <- forestedgeCovars[[c(18)]]
names(forestEdge) <- c("edge") # currently selecting 2022

# stacks for extracting values
covs = c(nonSeasonalCovars, distResamp, forestry, forestEdge)
covs_stack <- raster::stack(covs)
covs_stack <- as(covs_stack, "SpatRaster") # converting to spatraster

# quick field check
sort(names(covs_stack))
sort(colnames(selected_df[7:17]))


# ------------------------------------------------------------------------------
# DEFINING PREDICTION AND RASTERIZATION FUNCTIONS
# ------------------------------------------------------------------------------
pConstrain <- function(x, bounds){
  x <- ifelse(x < bounds[1], bounds[1], x)
  x <- ifelse(x > bounds[2], bounds[2], x)
  return(x)
}

# Function to extract covariate values at all cell locations in study area
gridCovs <- function(covs, st_area, puma_dat){
  
  # Get cell numbers for raster cells overlapping study area
  pred_cell_nums <- terra::cells(covs, vect(st_area))[,2]
  # Extract covariate values from across the covs raster for use in prediction surface
  pred_cov_temp <- terra::extract(covs, pred_cell_nums) 
  pred_xy <- terra::xyFromCell(covs, pred_cell_nums) %>% 
    as.data.frame() %>% 
    rename(utm_x = x, utm_y = y) 
  pred_covs <- cbind(pred_cov_temp, pred_xy) %>% 
    drop_na()
  
  # Rescale prediction covs based on range of values in rsf dataset
  covScale <- function(pc = pred_covs, dat = as.data.frame(puma_dat), cname){
    return((pc[,cname] - mean(dat[,cname], na.rm = TRUE))/sd(dat[,cname], na.rm = TRUE))
  }
  pred_covs_scale <- data.frame(agdist = covScale(cname = "ag_dist"),
                                agcov = covScale(cname = "all_ag_pcov_100m"),
                                devdist = covScale(cname = "developed_dist"),
                                devcov = covScale(cname = "developed_pcov_1km"),
                                activity = covScale(cname = "FPA_TH_activity"),
                                ripdist = covScale(cname = "riparian_dist"),
                                roaddist = covScale(cname = "road_dist"),
                                shrubcov = covScale(cname = "shrub_pcov_500m"),
                                slope = covScale(cname = "slope_1km"),
                                RAPtree = covScale(cname = "RAP_treecover"),
                                edge = covScale(cname = "edge")) %>% 
      mutate(agdist_sq = agdist^2,
             agcov_sq = agcov^2,
             devdist_sq = devdist^2,
             devcov_sq = devcov^2,
             roaddist_sq = roaddist^2,
             shrubcov_sq = shrubcov^2,
             slope_sq = slope^2,
             RAPtree_sq = RAPtree^2,
             edge_sq = edge^2,
             "devdist:slope" = devdist * slope,
             utm_x = pred_covs$utm_x,
             utm_y = pred_covs$utm_y
             )
  
  # Return scaled prediction covs
  return(pred_covs_scale)
}

# Function to predict from model equation at all cell locations in study area
gridPredict <- function(gcovs, mod, scale = "exp", constrain = NULL){
  coefs <- fixef(mod)[-c(1)] # removing intercept
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
gridRast <- function(dat, filename, scale = 30){
  # Convert prediction dataset to sf and then terra SpatVector
  pred_sp <- dat %>% 
    st_as_sf(coords = c("utm_x", "utm_y"), crs = su_crs) 
  
  # Create template raster
  r_temp <- rast()
  terra::ext(r_temp) <- terra::ext(pred_sp)
  terra::res(r_temp) <- scale
  terra::crs(r_temp) <- su_crs$wkt
  
  # Rasterize
  pred_sp <- as(pred_sp, "Spatial") # alternative coercing from sf via sp
  pred_sp <- vect(pred_sp)
  pred_rast <- terra::rasterize(pred_sp, r_temp, field = 'prediction', # original method
                                fun = "mean", filename = filename, overwrite = TRUE)
                                
  return(pred_rast)
}


# ------------------------------------------------------------------------------
# CREATE PREDICTION SURFACE
# ------------------------------------------------------------------------------
load(file = "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/top-final-age-sex-pooled-RSF.rda", verbose = TRUE) # load RSF
age_sex_pooled_covs = gridCovs(covs_stack, aoi, selected_df) 
age_sex_pooled_pred <- gridPredict(age_sex_pooled_covs, pooled_age_sex_models_final_top, scale = "exp") 
age_sex_pooled_rast <- gridRast(age_sex_pooled_pred, filename = "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/predictions/age-sex-pooled-exp-aoi-30.tif") 
