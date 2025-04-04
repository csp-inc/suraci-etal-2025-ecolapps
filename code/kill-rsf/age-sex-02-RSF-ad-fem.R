## ---------------------------
##
## Script name: age-sex-02-RSF-ad-fem.R
##
## Purpose: Fits RSF models for adult female pumas and then performs a model 
## selection routine
##
## Author: Mae Lacey
##
## Email contact: mae@csp-inc.org
##
## ---------------------------


library(ggeffects)
library(ggplot2)
library(raster)
library(purrr)
library(rgdal)
library(dplyr)
library(blme)

# Load model definitions
source('C:/Users/mae/Documents/github/pred-services/code/utils/rsf-age-sex-mod-select.R')
# Load kill site dataset
load("age-sex-class-RSFs/age_sex_RSF_covariate_extract.rda")
# Setting age/sex class of interest for modeling + file names
age_sex <- ad_fem
group <- "ad-fem"


# ------------------------------------------------------------------------------
# FIT AGE-SEX SUB-MODELS
# ------------------------------------------------------------------------------
age_sex_models <- list()
age_sex_models[[1]] <- ksModel(age_sex, human1)
age_sex_models[[2]] <- ksModel(age_sex, human2)
age_sex_models[[3]] <- ksModel(age_sex, human3)
age_sex_models[[4]] <- ksModel(age_sex, human4)
age_sex_models[[5]] <- ksModel(age_sex, human5)
age_sex_models[[6]] <- ksModel(age_sex, human6)
age_sex_models[[7]] <- ksModel(age_sex, human7) 
age_sex_models[[8]] <- ksModel(age_sex, human8)
age_sex_models[[9]] <- ksModel(age_sex, human9)
age_sex_models[[10]] <- ksModel(age_sex, human10)
age_sex_models[[11]] <- ksModel(age_sex, human11)
age_sex_models[[12]] <- ksModel(age_sex, human12)
age_sex_models[[13]] <- ksModel(age_sex, nh1)
age_sex_models[[14]] <- ksModel(age_sex, nh2)
age_sex_models[[15]] <- ksModel(age_sex, nh3)
age_sex_models[[16]] <- ksModel(age_sex, nh4)
age_sex_models[[17]] <- ksModel(age_sex, nh5)
age_sex_models[[18]] <- ksModel(age_sex, nh6)
age_sex_models[[19]] <- ksModel(age_sex, nh7)
age_sex_models[[20]] <- ksModel(age_sex, nh8)
age_sex_models[[21]] <- ksModel(age_sex, nh9)
age_sex_models[[22]] <- ksModel(age_sex, nh10)
age_sex_models[[23]] <- ksModel(age_sex, nh11)
age_sex_models[[24]] <- ksModel(age_sex, nh12)
age_sex_models[[25]] <- ksModel(age_sex, null_mod)

names(age_sex_models) <- c('human1','human2','human3','human4','human5', 'human6',
                           'human7','human8','human9','human10','human11','human12',
                           'nh1','nh2','nh3','nh4','nh5','nh6','nh7','nh8','nh9',
                           'nh10','nh11','nh12','null_mod')

# Save outputs
save(list = c("age_sex","age_sex_models"), 
     file = paste0("age-sex-class-RSFs/model-outputs/prelim-RSFs-", group, ".rda"))

# Create model selection table, exclude models w/ convergence issues
mod_select <- purrr::map(age_sex_models, function(x){!is.na(AIC(x)) & !is.nan(AIC(x))}) %>% unlist()
age_sex_models_good <- age_sex_models[mod_select] 
age_sex_models_tab <- bbmle::ICtab(age_sex_models_good, mnames = names(age_sex_models_good), 
                                  type=c("AIC","BIC","AICc","qAIC","qAICc"),
                                  weights = FALSE, delta = TRUE, base = TRUE,
                                  logLik=FALSE, sort = TRUE,
                                  nobs=NULL, dispersion = 1, k = 2)
age_sex_models_tab <- as.data.frame(age_sex_models_tab)
# save top model
age_sex_models_top_human <- age_sex_models$human8
age_sex_models_top_nh <- age_sex_models$nh3
save(list = c('age_sex_models_tab', 'age_sex_models_top_human', 'age_sex_models_top_nh'), 
     file = paste0("age-sex-class-RSFs/model-outputs/top-prelim-RSFs-", group, ".rda"))
write.csv(age_sex_models_tab, 
          paste0("age-sex-class-RSFs/model-outputs/prelim-RSFs-AICoutputs-", group, ".csv"))


# ------------------------------------------------------------------------------
# FIT TOP AGE-SEX MODELS + COMBOS
# ------------------------------------------------------------------------------
age_sex_models_final <- list()
age_sex_models_final[[1]] <- ksModel(age_sex, human8)
age_sex_models_final[[2]] <- ksModel(age_sex, nh3)
age_sex_models_final[[3]] <- ksModel(age_sex, af_combo1)
age_sex_models_final[[4]] <- ksModel(age_sex, af_combo2)
age_sex_models_final[[5]] <- ksModel(age_sex, null_mod)

names(age_sex_models_final) <- c('human8','nh3','af_combo1', 'af_combo2', 'null_mod')

# Save outputs
save(list = c("age_sex","age_sex_models_final"), 
     file = paste0("age-sex-class-RSFs/model-outputs/final-RSFs-", group, ".rda"))

# Create model selection table, exclude models w/ convergence issues
mod_select_final <- purrr::map(age_sex_models_final, function(x){!is.na(AIC(x)) & !is.nan(AIC(x))}) %>% unlist()
age_sex_models_final_good <- age_sex_models_final[mod_select_final] 
age_sex_models_final_tab <- bbmle::ICtab(age_sex_models_final_good, mnames = names(age_sex_models_final_good), 
                                  type=c("AIC","BIC","AICc","qAIC","qAICc"),
                                  weights = FALSE, delta = TRUE, base = TRUE,
                                  logLik=FALSE, sort = TRUE,
                                  nobs=NULL, dispersion = 1, k = 2)
age_sex_models_final_tab <- as.data.frame(age_sex_models_final_tab)
# save top model
age_sex_models_final_top <- age_sex_models_final$af_combo2
save(list = c('age_sex_models_final_tab', 'age_sex_models_final_top'), 
     file = paste0("age-sex-class-RSFs/model-outputs/top-final-RSF-", group, ".rda"))
write.csv(age_sex_models_final_tab, 
          paste0("age-sex-class-RSFs/model-outputs/final-RSFs-AICoutputs-", group, ".csv"))