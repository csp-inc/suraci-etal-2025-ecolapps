## ---------------------------
##
## Script name: age-sex-pooled-01-enviro-covar-extract-SOE-select.R
##
## Purpose: Cleans kill site data, extracts environmental covariates at all kill
## site locations, then performs a scale of effect analysis to select the relevant 
## scale at which each covariate is processed. A final dataframe retaining covariates
## processed at the scales selected is then prepared and exported for modeling.
##
## Author: Mae Lacey
##
## Email contact: mae@csp-inc.org
##
## ---------------------------


###################################################
## Clean up kill site data and remove duplicates ##
###################################################

library(tidyverse)
library(corrplot)
library(stringr)
library(ggplot2)
library(ggpubr)
library(raster)
library(repmod)
library(dplyr)
library(lme4)
library(sf)
library(AICcmodavg)
library(bbmle)


# ------------------------------------------------------------------------------
# READ IN & PROCESS POINT DATA
# ------------------------------------------------------------------------------
# we'll be running the scale of effect analysis on the full dataset
# load in cleaned kill site data w/ 10 random "available" points for each kill site
ks_bg_pt <- read.csv("age-sex-class-RSFs/analysis-updates_kill-used-available-age-sex-dataset-20231212.csv") 
# join prey type from original dataset back to clean dataset for prey-specific analyses
prey <- read.csv("kill-sites-20221031.csv") %>% 
  dplyr::select(Report_Id, Prey_Species, FinalPrey_LMLadded)
ks_bg <- left_join(ks_bg_pt, prey, by = "Report_Id")
# convert to spatial data frame
coordinates(ks_bg) <- ~UTM_X+UTM_Y
crs(ks_bg) = "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs"


# ------------------------------------------------------------------------------
# READ IN & PROCESS RASTERS
# ------------------------------------------------------------------------------
# Load all rasters for processing
# (note that we're eliminating all seasonal covs (i.e., climate covs) for age-sex analysis)
nsResamp <- stack("D:/Panthera/RawData/final-covars/nsResamp.grd") # version of 'raw-covs/puma-kill-site-covs-non-seasonal-90m.tif' resampled to 30m
distResamp <- stack("D:/Panthera/RawData/final-covars/distResamp.grd") # version of 'raw-covs/puma-kill-site-covs-distvars-270m.tif' resampled to 30m
forestCovars <- stack("raw-covs/puma-kill-site-covs-forestvars-30m.tif")
forestedgeCovars <- stack("raw-covs/puma-kill-site-covs-forestedgevars-30m.tif")

# Select all relevant layers
# non-seasonal covariates - dropping VRM, NLCD forest % cover, impervious, human pop den, and nightlight
nsResamp <- dropLayer(nsResamp, c(4,8,12,13,16,19,22,23,24,25,32,33,34,35,36,37,38,39))
# forestry covariates - only keep smoothed covariates (nothing to remove from forestedgeCovars or distCovars)
forestry <- forestCovars[[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,25)]]

# Confirm that all extents match
extent(nsResamp)
extent(distResamp)
extent(forestedgeCovars)
extent(forestry)


# ------------------------------------------------------------------------------
# CREATE FINAL COV STACK & EXTRACT VALUES AT POINTS
# ------------------------------------------------------------------------------
# stack together
all_covs <- stack(nsResamp, distResamp, forestedgeCovars, forestry)
names(all_covs) # confirm names

# extract raster values at points
cov_extract <- raster::extract(all_covs, ks_bg)
cov_extract_df <- as.data.frame(cbind(ks_bg, cov_extract))
save(cov_extract_df, file = "age-sex-class-RSFs/age-sex-pooled-RSF/age_sex_pooled_RSF_covariate_extract_for_SOE.rda")


# ------------------------------------------------------------------------------
# SCALE OF EFFECT ANALYSIS:
# ------------------------------------------------------------------------------
# for all covariates processed at three radii, run a single-variable model to select the most appropriate radius

# prep & set up empty lists to add models to iteratively
covars <- colnames(cov_extract_df[10:74]) # select only covariate fields         
sorted_covs <- covars[order(unlist(covars),decreasing=FALSE)]
mod_list <- list()
mod_names <- list()

# run single-variable models
for (i in 1:length(sorted_covs)){                                                       
  print(sorted_covs[i])
  mod_eq <- paste0("kill ~ ", sorted_covs[i], " + (1 | Cougar_ID)") %>% 
    as.formula()
  print(mod_eq)
  mod_temp <- glmer(mod_eq, data = cov_extract_df, family = "binomial")
  mod_list[[i]] <- mod_temp
  mod_names[[i]] <- sorted_covs[i]
}
names(mod_list) <- mod_names

# create model selection table, exclude models w/ convergence issues
mod_select <- purrr::map(mod_list, function(x){!is.na(AIC(x)) & !is.nan(AIC(x))}) %>% unlist()
mod_list_good <- mod_list[mod_select] 
mod_list_tab <- bbmle::ICtab(mod_list_good, mnames = names(mod_list_good), 
                                  type=c("AIC","BIC","AICc","qAIC","qAICc"),
                                  weights = FALSE, delta = TRUE, base = TRUE,
                                  logLik=FALSE, sort = TRUE,
                                  nobs=NULL, dispersion = 1, k = 2)
mod_list_tab <- as.data.frame(mod_list_tab)
mod_list_tab$mods <- rownames(mod_list_tab) # add rownames to column
save(mod_list, mod_list_tab, file = "age-sex-class-RSFs/age-sex-pooled-RSF/age_sex_pooled_RSF_SOE_models.rda")

# retain only one version of each covariate, based on lowest AIC from SOE analysis
selected_df <- cov_extract_df %>%
  dplyr::select( 
    Cougar_ID,
    Report_Id,
    age_sex,
    year,
    kill,
    FinalPrey_LMLadded,
    ag_dist,
    all_ag_pcov_100m,
    aspect_100m, 
    developed_dist,
    developed_pcov_1km,
    edge_2017_pcov_100m,
    edge_2018_pcov_100m,
    edge_2019_pcov_100m,
    edge_2020_pcov_100m,
    edge_2021_pcov_100m,
    edge_2022_pcov_100m,
    elevation_100m, 
    FPA_TH_activity,
    FPA_TH_invert_forestmaskpostedge_disturbed_edge_pcov_1km,
    grass_pcov_1km, 
    riparian_dist,
    road_dist, 
    shrub_pcov_500m, 
    slope_1km,
    TRE_2017_pcov_1km,
    TRE_2018_pcov_1km,
    TRE_2019_pcov_1km,
    TRE_2020_pcov_1km,
    TRE_2021_pcov_1km,
    TRE_2022_pcov_1km
  )

# rename one covariate w/ long name
colnames(selected_df)[colnames(selected_df) == "FPA_TH_invert_forestmaskpostedge_disturbed_edge_pcov_1km"] ="disturbed_pcov_1km"


# ------------------------------------------------------------------------------
# CREATE ONE COLUMN FOR EACH ANNUAL DATASET
# ------------------------------------------------------------------------------
# year match covariates to points
selected_df <- selected_df %>% mutate(
  RAP_treecover = case_when(year == 2017 ~ TRE_2017_pcov_1km,
                            year == 2018 ~ TRE_2018_pcov_1km,
                            year == 2019 ~ TRE_2019_pcov_1km,
                            year == 2020 ~ TRE_2020_pcov_1km,
                            year == 2021 ~ TRE_2021_pcov_1km,
                            year == 2022 ~ TRE_2022_pcov_1km),
  edge = case_when(year == 2017 ~ edge_2017_pcov_100m,
                   year == 2018 ~ edge_2018_pcov_100m,
                   year == 2019 ~ edge_2019_pcov_100m,
                   year == 2020 ~ edge_2020_pcov_100m,
                   year == 2021 ~ edge_2021_pcov_100m,
                   year == 2022 ~ edge_2022_pcov_100m)) %>%
  dplyr::select(-c(TRE_2017_pcov_1km, 
                   TRE_2018_pcov_1km,
                   TRE_2019_pcov_1km,
                   TRE_2020_pcov_1km,
                   TRE_2021_pcov_1km,
                   TRE_2022_pcov_1km,
                   edge_2017_pcov_100m,
                   edge_2018_pcov_100m,
                   edge_2019_pcov_100m,
                   edge_2020_pcov_100m,
                   edge_2021_pcov_100m,
                   edge_2022_pcov_100m
                   ))

save(selected_df, file = "age-sex-class-RSFs/age-sex-pooled-RSF/age_sex_pooled_RSF_covariate_extract_all.rda")


# ---------------------------------------------------------------------------
# CHECK FOR CORRELATIONS IN FINAL DATASET 
# ---------------------------------------------------------------------------
# select final covariates based on AIC of single-variable models above
corr_df <- selected_df %>%
  dplyr::select(-c(Cougar_ID, 
                   Report_Id,
                   age_sex,
                   year,
                   kill,
                   FinalPrey_LMLadded))

# correlation
final_corr <- as.data.frame(cor(corr_df, method = c("pearson"), use = "complete.obs"))
#final_cplot <- corrplot(final_corr, type = "upper", order = "hclust", 
#                               tl.col = "black", tl.srt = 45)
write.csv(final_corr, "age-sex-class-RSFs/age-sex-pooled-RSF/age_sex_pooled_RSF_final_corr.csv") 
save(selected_df, mod_list_tab, final_corr, file = "age-sex-class-RSFs/age-sex-pooled-RSF/age_sex_pooled_RSF_covExtract_SOE_correlations.rda")


# ---------------------------------------------------------------------------
# CREATE FINAL DATAFRAME FOR ANALYSIS
# ---------------------------------------------------------------------------
all_age_sex <- selected_df %>%
  mutate(agdist = c(scale(ag_dist)),
         agdist_sq = agdist^2,
         agcov = c(scale(all_ag_pcov_100m)),
         agcov_sq = agcov^2,
         aspect = c(scale(aspect_100m)),
         devdist = c(scale(developed_dist)),
         devdist_sq = devdist^2,
         devcov = c(scale(developed_pcov_1km)),
         devcov_sq = devcov^2,
         activity = c(scale(FPA_TH_activity)),
         ripdist = c(scale(riparian_dist)),
         roaddist = c(scale(road_dist)),
         roaddist_sq = roaddist^2,        
         shrubcov = c(scale(shrub_pcov_500m)),
         shrubcov_sq = shrubcov^2,
         slope = c(scale(slope_1km)),
         slope_sq = slope^2,
         RAPtree = c(scale(RAP_treecover)),
         RAPtree_sq = RAPtree^2,
         edge = c(scale(edge)),
         edge_sq = edge^2) %>%
  dplyr::select(Cougar_ID, Report_Id, age_sex, year, kill, agdist, agdist_sq,
                agcov, agcov_sq, aspect, devdist, devdist_sq, devcov, devcov_sq,
                activity, ripdist, roaddist, roaddist_sq, shrubcov, shrubcov_sq,   
                slope, slope_sq, RAPtree, RAPtree_sq, edge, edge_sq) 

# save age-sex classes:
save(list = c("all_age_sex"), 
     file = "age-sex-class-RSFs/age-sex-pooled-RSF/age_sex_pooled_RSF_covariate_extract.rda")


# ------------------------------------------------------------------------------
# GENERATE PLOTS OF COVARIATE DISTRIBUTIONS FOR INSPECTION
# ------------------------------------------------------------------------------
plot1 <- ggplot(corr_df, aes(x=aspect_100m)) + geom_histogram(binwidth=10)
plot2 <- ggplot(corr_df, aes(x=elevation_100m)) + geom_histogram(binwidth=10)
plot3 <- ggplot(corr_df, aes(x=grass_pcov_1km)) + geom_histogram(binwidth=0.01)
plot4 <- ggplot(corr_df, aes(x=slope_1km)) + geom_histogram(binwidth=1)
plot5 <- ggplot(corr_df, aes(x=shrub_pcov_500m)) + geom_histogram(binwidth=0.01)
plot6 <- ggplot(corr_df, aes(x=riparian_dist)) + geom_histogram(binwidth=10)
plot7 <- ggplot(corr_df, aes(x=road_dist)) + geom_histogram(binwidth=50)
plot8 <- ggplot(corr_df, aes(x=ag_dist)) + geom_histogram(binwidth=10)
plot9 <- ggplot(corr_df, aes(x=all_ag_pcov_100m)) + geom_histogram(binwidth=0.01)
plot10 <- ggplot(corr_df, aes(x=developed_pcov_1km)) + geom_histogram(binwidth=0.01)
plot11 <- ggplot(corr_df, aes(x=developed_dist)) + geom_histogram(binwidth=10)
plot12 <- ggplot(corr_df, aes(x=FPA_TH_activity)) + geom_histogram(binwidth=1)
plot13 <- ggplot(corr_df, aes(x=disturbed_pcov_1km)) + geom_histogram(binwidth=0.01)
plot14 <- ggplot(corr_df, aes(x=RAP_treecover)) + geom_histogram(binwidth=1)
plot15 <- ggplot(corr_df, aes(x=edge)) + geom_histogram(binwidth=0.1)


plots <- ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, 
                   plot11, plot12, plot13, plot14, plot15,
                   ncol = 4, nrow = 5)
plots
# then export
ggsave(plots, file="age-sex-class-RSFs/age-sex-pooled-RSF/age_sex_pooled_RSF_final_cov_plots.png",
       width=10, height=7)