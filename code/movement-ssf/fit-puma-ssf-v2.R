library(tidyverse)
library(sf)
library(terra)
library(amt)
library(lubridate)
library(glmmTMB)
library(MuMIn)
library(AICcmodavg)

# Load model definitions
source('code/utils/ssf-mod-lkp-v2.R')

# ---------- Compile Analysis Datasets -------------------
load(file = "data/puma-ssf-dataset-20230701.rda", verbose = TRUE)
ssf_dat <- puma_ssf_dat %>% 
  dplyr::mutate(case = as.numeric(case_))

#------------- Analysis Functions ------------------
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

modFit <- function(data, model){
  mod <- glmmTMB(model, data = data, family = poisson, doFit = FALSE,
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
  mod$parameters$theta[1] <- log(1e3)
  ntheta = length(mod$parameters$theta)
  mod$mapArg <- list(theta=factor(c(NA, rep(1,ntheta-1))))
  mod_fit <- glmmTMB:::fitTMB(mod)
  return(mod_fit)
}


#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# STAGE 1 - SUB-MODEL SELECTION

#---------------------------------------------
#------------- ADULT FEMALE ------------------
#---------------------------------------------
# Prep dataset
ad_fem <- ssf_dat %>% filter(age_class == "adult", sex == "Female") %>% datPrep()

# MODEL SELECTION - ROUND 1
# Fit models -HUMAN
ad_fem_hmods <- list()
ad_fem_hmods[[1]] <- modFit(ad_fem, human1)
ad_fem_hmods[[2]] <- modFit(ad_fem, human2)
ad_fem_hmods[[3]] <- modFit(ad_fem, human3)
ad_fem_hmods[[4]] <- modFit(ad_fem, human4)
ad_fem_hmods[[5]] <- modFit(ad_fem, human5)
ad_fem_hmods[[6]] <- modFit(ad_fem, human6)
ad_fem_hmods[[7]] <- modFit(ad_fem, human7)
ad_fem_hmods[[8]] <- modFit(ad_fem, human8)
ad_fem_hmods[[9]] <- modFit(ad_fem, human9)
ad_fem_hmods[[10]] <- modFit(ad_fem, human10)

names(ad_fem_hmods) <- c('human1','human2','human3','human4',
                         'human5','human6','human7','human8',
                         'human9','human10')

# Fit models -NON-HUMAN
ad_fem_nhmods <- list()
ad_fem_nhmods[[1]] <- modFit(ad_fem, nh1)
ad_fem_nhmods[[2]] <- modFit(ad_fem, nh2)
ad_fem_nhmods[[3]] <- modFit(ad_fem, nh3)
ad_fem_nhmods[[4]] <- modFit(ad_fem, nh4)
ad_fem_nhmods[[5]] <- modFit(ad_fem, nh5)
ad_fem_nhmods[[6]] <- modFit(ad_fem, nh6)
ad_fem_nhmods[[7]] <- modFit(ad_fem, nh7)
ad_fem_nhmods[[8]] <- modFit(ad_fem, nh8)

names(ad_fem_nhmods) <- c('nh1','nh2','nh3','nh4',
                          'nh5','nh6','nh7','nh8')

# Save everything
save(list = c("ad_fem","ad_fem_hmods","ad_fem_nhmods"), file = "output/ad-fem-ssf-mods-v2.rda")

#---------------------------------------------
#------------- ADULT MALE ------------------
#---------------------------------------------
# Prep dataset
ad_male <- ssf_dat %>% filter(age_class == "adult", sex == "Male") %>% datPrep()

# MODEL SELECTION - ROUND 1
# Fit models -HUMAN
ad_male_hmods <- list()
ad_male_hmods[[1]] <- modFit(ad_male, human1)
ad_male_hmods[[2]] <- modFit(ad_male, human2)
ad_male_hmods[[3]] <- modFit(ad_male, human3)
ad_male_hmods[[4]] <- modFit(ad_male, human4)
ad_male_hmods[[5]] <- modFit(ad_male, human5)
ad_male_hmods[[6]] <- modFit(ad_male, human6)
ad_male_hmods[[7]] <- modFit(ad_male, human7)
ad_male_hmods[[8]] <- modFit(ad_male, human8)
ad_male_hmods[[9]] <- modFit(ad_male, human9)
ad_male_hmods[[10]] <- modFit(ad_male, human10)

names(ad_male_hmods) <- c('human1','human2','human3','human4',
                         'human5','human6','human7','human8',
                         'human9','human10')

# Fit models -NON-HUMAN
ad_male_nhmods <- list()
ad_male_nhmods[[1]] <- modFit(ad_male, nh1)
ad_male_nhmods[[2]] <- modFit(ad_male, nh2)
ad_male_nhmods[[3]] <- modFit(ad_male, nh3)
ad_male_nhmods[[4]] <- modFit(ad_male, nh4)
ad_male_nhmods[[5]] <- modFit(ad_male, nh5)
# ad_male_nhmods[[6]] <- modFit(ad_male, nh6)
ad_male_nhmods[[6]] <- modFit(ad_male, nh7)
ad_male_nhmods[[7]] <- modFit(ad_male, nh8)

names(ad_male_nhmods) <- c('nh1','nh2','nh3','nh4',
                           'nh5','nh7','nh8')

# Save everything
save(list = c("ad_male","ad_male_hmods","ad_male_nhmods"), file = "output/ad-male-ssf-mods-v2.rda")


#---------------------------------------------
#------------- DISPERSER FEMALE --------------
#---------------------------------------------
# Prep dataset
disp_fem <- ssf_dat %>% filter(age_class == "disperser", sex == "Female") %>% datPrep()

# MODEL SELECTION - ROUND 1
# Fit models - HUMAN
disp_fem_hmods <- list()
disp_fem_hmods[[1]] <- modFit(disp_fem, human1)
disp_fem_hmods[[2]] <- modFit(disp_fem, human2)
disp_fem_hmods[[3]] <- modFit(disp_fem, human3)
disp_fem_hmods[[4]] <- modFit(disp_fem, human4)
disp_fem_hmods[[5]] <- modFit(disp_fem, human5)
disp_fem_hmods[[6]] <- modFit(disp_fem, human6)
disp_fem_hmods[[7]] <- modFit(disp_fem, human7)
disp_fem_hmods[[8]] <- modFit(disp_fem, human8)
disp_fem_hmods[[9]] <- modFit(disp_fem, human9)
disp_fem_hmods[[10]] <- modFit(disp_fem, human10)

names(disp_fem_hmods) <- c('human1','human2','human3','human4',
                         'human5','human6','human7','human8',
                         'human9','human10')

# Fit models -NON-HUMAN
disp_fem_nhmods <- list()
disp_fem_nhmods[[1]] <- modFit(disp_fem, nh1)
disp_fem_nhmods[[2]] <- modFit(disp_fem, nh2)
disp_fem_nhmods[[3]] <- modFit(disp_fem, nh3)
disp_fem_nhmods[[4]] <- modFit(disp_fem, nh4)
disp_fem_nhmods[[5]] <- modFit(disp_fem, nh5)
disp_fem_nhmods[[6]] <- modFit(disp_fem, nh6)
disp_fem_nhmods[[7]] <- modFit(disp_fem, nh7)
disp_fem_nhmods[[8]] <- modFit(disp_fem, nh8)

names(disp_fem_nhmods) <- c('nh1','nh2','nh3',
                          'nh4','nh5','nh6','nh7','nh8')

# Save everything
save(list = c("disp_fem","disp_fem_hmods","disp_fem_nhmods"), file = "output/disp-fem-ssf-mods-v2.rda")

#---------------------------------------------
#------------- DISPERSER MALE ----------------
#---------------------------------------------
# Prep dataset
disp_male <- ssf_dat %>% filter(age_class == "disperser", sex == "Male") %>% datPrep()

# MODEL SELECTION - ROUND 1
# Fit models -HUMAN
disp_male_hmods <- list()
disp_male_hmods[[1]] <- modFit(disp_male, human1)
disp_male_hmods[[2]] <- modFit(disp_male, human2)
disp_male_hmods[[3]] <- modFit(disp_male, human3)
disp_male_hmods[[4]] <- modFit(disp_male, human4)
disp_male_hmods[[5]] <- modFit(disp_male, human5)
disp_male_hmods[[6]] <- modFit(disp_male, human6)
disp_male_hmods[[7]] <- modFit(disp_male, human7)
disp_male_hmods[[8]] <- modFit(disp_male, human8)
disp_male_hmods[[9]] <- modFit(disp_male, human9)
disp_male_hmods[[10]] <- modFit(disp_male, human10)

names(disp_male_hmods) <- c('human1','human2','human3','human4',
                          'human5','human6','human7','human8',
                          'human9','human10')

# Fit models -NON-HUMAN
disp_male_nhmods <- list()
disp_male_nhmods[[1]] <- modFit(disp_male, nh1)
disp_male_nhmods[[2]] <- modFit(disp_male, nh2)
disp_male_nhmods[[3]] <- modFit(disp_male, nh3)
disp_male_nhmods[[4]] <- modFit(disp_male, nh4)
disp_male_nhmods[[5]] <- modFit(disp_male, nh5)
disp_male_nhmods[[6]] <- modFit(disp_male, nh6)
disp_male_nhmods[[7]] <- modFit(disp_male, nh7)
disp_male_nhmods[[8]] <- modFit(disp_male, nh8)

names(disp_male_nhmods) <- c('nh1','nh2','nh3',
                           'nh4','nh5','nh6',
                           'nh7','nh8')

# Exclude mods with convergence issues
g <- map(disp_male_nhmods, function(x){!is.na(AIC(x)) & !is.nan(AIC(x))}) %>% unlist()
disp_male_nhmods <- disp_male_nhmods[g]

# Save everything
save(list = c("disp_male","disp_male_hmods","disp_male_nhmods"), file = "output/disp-male-ssf-mods-v2.rda")


#---------------------------------------------
#--------------- ALL CATS --------------------
#---------------------------------------------
# Prep dataset
all_cats <- ssf_dat %>% datPrep()

# MODEL SELECTION - ROUND 1
# Fit models -HUMAN
all_cats_hmods <- list()
all_cats_hmods[[1]] <- modFit(all_cats, human1)
all_cats_hmods[[2]] <- modFit(all_cats, human2)
all_cats_hmods[[3]] <- modFit(all_cats, human3)
all_cats_hmods[[4]] <- modFit(all_cats, human4)
all_cats_hmods[[5]] <- modFit(all_cats, human5)
all_cats_hmods[[6]] <- modFit(all_cats, human6)
all_cats_hmods[[7]] <- modFit(all_cats, human7)
all_cats_hmods[[8]] <- modFit(all_cats, human8)
all_cats_hmods[[9]] <- modFit(all_cats, human9)
all_cats_hmods[[10]] <- modFit(all_cats, human10)

names(all_cats_hmods) <- c('human1','human2','human3','human4',
                         'human5','human6','human7','human8',
                         'human9','human10')

# Fit models -NON-HUMAN
all_cats_nhmods <- list()
all_cats_nhmods[[1]] <- modFit(all_cats, nh1)
all_cats_nhmods[[2]] <- modFit(all_cats, nh2)
all_cats_nhmods[[3]] <- modFit(all_cats, nh3)
all_cats_nhmods[[4]] <- modFit(all_cats, nh4)
all_cats_nhmods[[5]] <- modFit(all_cats, nh5)
all_cats_nhmods[[6]] <- modFit(all_cats, nh6)
all_cats_nhmods[[7]] <- modFit(all_cats, nh7)
all_cats_nhmods[[8]] <- modFit(all_cats, nh8)

names(all_cats_nhmods) <- c('nh1','nh2','nh3','nh4',
                          'nh5','nh6','nh7','nh8')

# Save everything
save(list = c("all_cats","all_cats_hmods","all_cats_nhmods"), file = "output/all-cats-ssf-mods.rda")

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# STAGE 2 - COMBO MODEL SELECTION

#---------------------------------------------
#------------- ADULT FEMALE ------------------
#---------------------------------------------
load(file = "output/ad-fem-ssf-mods-v2.rda", verbose = T)

# Identify top human and non-human sub models
af_tab_human<-aictab(ad_fem_hmods)
af_tab_nh<-aictab(ad_fem_nhmods)

# Fit top sub models and combos
ad_fem_combo <- list()
ad_fem_combo[[1]] <- ad_fem_hmods$human5
ad_fem_combo[[2]] <- ad_fem_nhmods$nh2
ad_fem_combo[[3]] <- modFit(ad_fem, af_combo1)
ad_fem_combo[[4]] <- modFit(ad_fem, af_combo2)
ad_fem_combo[[5]] <- modFit(ad_fem, af_combo3)

names(ad_fem_combo) <- c('topHuman','topNonhuman','combo1','combo2','combo3')

# Identify top mod overall
af_tab_combo<-aictab(ad_fem_combo)
ad_fem_top <- ad_fem_combo$combo2

#save
save(list = c('af_tab_human', 'af_tab_nh', 'af_tab_combo', 'ad_fem_top'), file = "output/ad-fem-ssf-top-mod-v2.rda")

#---------------------------------------------
#------------- ADULT MALE ------------------
#---------------------------------------------
load(file = "output/ad-male-ssf-mods-v2.rda", verbose = T)

# Identify top human and non-human sub models
am_tab_human <- aictab(ad_male_hmods)
am_tab_nh <- aictab(ad_male_nhmods)

# Fit top sub models and combos
ad_male_combo <- list()
ad_male_combo[[1]] <- ad_male_hmods$human5
ad_male_combo[[2]] <- ad_male_nhmods$nh2
ad_male_combo[[3]] <- modFit(ad_male, am_combo1)
ad_male_combo[[4]] <- modFit(ad_male, am_combo2)

names(ad_male_combo) <- c('topHuman','topNonhuman','combo1','combo2')

# Identify top mod overall
am_tab_combo<-aictab(ad_male_combo)
ad_male_top <- ad_male_combo$combo2

# Save
save(list = c('am_tab_human', 'am_tab_nh', 'am_tab_combo', 'ad_male_top'), file = "output/ad-male-ssf-top-mod-v2.rda")

#---------------------------------------------
#------------- DISPERSER FEMALE ------------------
#---------------------------------------------
load(file = "output/disp-fem-ssf-mods-v2.rda", verbose = T)

# Identify top human and non-human sub models
df_tab_human<-aictab(disp_fem_hmods)
df_tab_nh<-aictab(disp_fem_nhmods)

# Fit top sub models and combos
disp_fem_combo <- list()
disp_fem_combo[[1]] <- disp_fem_hmods$human7
disp_fem_combo[[2]] <- disp_fem_nhmods$nh2
disp_fem_combo[[3]] <- modFit(disp_fem, df_combo1)
disp_fem_combo[[4]] <- modFit(disp_fem, df_combo2)
disp_fem_combo[[5]] <- modFit(disp_fem, df_combo3)

names(disp_fem_combo) <- c('topHuman','topNonhuman','combo1','combo2','combo3')

# Identify top mod overall
df_tab_combo<-aictab(disp_fem_combo)
disp_fem_top <- disp_fem_combo$combo2

# Save
save(list = c('df_tab_human', 'df_tab_nh', 'df_tab_combo', 'disp_fem_top'), file = "output/disp-fem-ssf-top-mod-v2.rda")

#---------------------------------------------
#------------- DISPERSER MALE ------------------
#---------------------------------------------
load(file = "output/disp-male-ssf-mods-v2.rda", verbose = T)

# Identify top human and non-human sub models
dm_tab_human <- aictab(disp_male_hmods)
dm_tab_nh <- aictab(disp_male_nhmods)

# Fit top sub models and combos
disp_male_combo <- list()
disp_male_combo[[1]] <- disp_male_hmods$human5
disp_male_combo[[2]] <- disp_male_nhmods$nh2
disp_male_combo[[3]] <- modFit(disp_male, dm_combo1)
disp_male_combo[[4]] <- modFit(disp_male, dm_combo2)

names(disp_male_combo) <- c('topHuman','topNonhuman','combo1','combo2')

#---------------------------------------------
#--------------- ALL CATS --------------------
#---------------------------------------------
load(file = "output/all-cats-ssf-mods.rda", verbose = T)

# Identify top human and non-human sub models
ac_tab_human<-aictab(all_cats_hmods[-4])
ac_tab_nh<-aictab(all_cats_nhmods[-7])

# Fit top sub models and combos
all_cats_combo <- list()
all_cats_combo[[1]] <- all_cats_hmods$human5
all_cats_combo[[2]] <- all_cats_nhmods$nh2
all_cats_combo[[3]] <- modFit(all_cats, ac_combo1)
all_cats_combo[[4]] <- modFit(all_cats, ac_combo2)
# all_cats_combo[[5]] <- modFit(all_cats, ac_combo3)

names(all_cats_combo) <- c('topHuman','topNonhuman','combo1','combo2')

# Identify top mod overall
ac_tab_combo<-aictab(all_cats_combo)
all_cats_top <- all_cats_combo$combo1

#save
save(list = c('ac_tab_human', 'ac_tab_nh', 'ac_tab_combo', 'all_cats_top'), file = "output/all-cats-ssf-top-mod.rda")
