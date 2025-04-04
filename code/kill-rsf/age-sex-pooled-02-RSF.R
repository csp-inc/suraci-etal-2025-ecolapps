## ---------------------------
##
## Script name: age-sex-pooled-02-RSF.R
##
## Purpose: Fits pooled RSF models and then performs a model selection routine
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
load("age-sex-class-RSFs/age-sex-pooled-RSF/age_sex_pooled_RSF_covariate_extract.rda")
# Setting age/sex class of interest for modeling + file names
age_sex <- all_age_sex # here, we're using the full dataset (as opposed to specific age-sex groups)


# ------------------------------------------------------------------------------
# FIT AGE-SEX SUB-MODELS
# ------------------------------------------------------------------------------
pooled_age_sex_models <- list()
pooled_age_sex_models[[1]] <- ksModel(age_sex, human1)
pooled_age_sex_models[[2]] <- ksModel(age_sex, human2)
pooled_age_sex_models[[3]] <- ksModel(age_sex, human3)
pooled_age_sex_models[[4]] <- ksModel(age_sex, human4)
pooled_age_sex_models[[5]] <- ksModel(age_sex, human5)
pooled_age_sex_models[[6]] <- ksModel(age_sex, human6)
pooled_age_sex_models[[7]] <- ksModel(age_sex, human7) 
pooled_age_sex_models[[8]] <- ksModel(age_sex, human8)
pooled_age_sex_models[[9]] <- ksModel(age_sex, human9)
pooled_age_sex_models[[10]] <- ksModel(age_sex, human10)
pooled_age_sex_models[[11]] <- ksModel(age_sex, human11)
pooled_age_sex_models[[12]] <- ksModel(age_sex, human12)
pooled_age_sex_models[[13]] <- ksModel(age_sex, nh1)
pooled_age_sex_models[[14]] <- ksModel(age_sex, nh2)
pooled_age_sex_models[[15]] <- ksModel(age_sex, nh3)
pooled_age_sex_models[[16]] <- ksModel(age_sex, nh4)
pooled_age_sex_models[[17]] <- ksModel(age_sex, nh5)
pooled_age_sex_models[[18]] <- ksModel(age_sex, nh6)
pooled_age_sex_models[[19]] <- ksModel(age_sex, nh7)
pooled_age_sex_models[[20]] <- ksModel(age_sex, nh8)
pooled_age_sex_models[[21]] <- ksModel(age_sex, nh9)
pooled_age_sex_models[[22]] <- ksModel(age_sex, nh10)
pooled_age_sex_models[[23]] <- ksModel(age_sex, nh11)
pooled_age_sex_models[[24]] <- ksModel(age_sex, nh12)
pooled_age_sex_models[[25]] <- ksModel(age_sex, null_mod)

names(pooled_age_sex_models) <- c('human1','human2','human3','human4','human5', 'human6',
                           'human7','human8','human9','human10','human11','human12',
                           'nh1','nh2','nh3','nh4','nh5','nh6','nh7','nh8','nh9',
                           'nh10','nh11','nh12','null_mod')

# Save outputs
save(list = c("age_sex","pooled_age_sex_models"), 
     file = "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/prelim-age-sex-pooled-RSFs.rda")

# Create model selection table, exclude models w/ convergence issues
mod_select <- purrr::map(pooled_age_sex_models, function(x){!is.na(AIC(x)) & !is.nan(AIC(x))}) %>% unlist()
pooled_age_sex_models_good <- pooled_age_sex_models[mod_select] 
pooled_age_sex_models_tab <- bbmle::ICtab(pooled_age_sex_models_good, mnames = names(pooled_age_sex_models_good), 
                                   type=c("AIC","BIC","AICc","qAIC","qAICc"),
                                   weights = FALSE, delta = TRUE, base = TRUE,
                                   logLik=FALSE, sort = TRUE,
                                   nobs=NULL, dispersion = 1, k = 2)
pooled_age_sex_models_tab <- as.data.frame(pooled_age_sex_models_tab)
# save top model
pooled_age_sex_models_top_human <- pooled_age_sex_models$human6
pooled_age_sex_models_top_nh <- pooled_age_sex_models$nh11
save(list = c('pooled_age_sex_models_tab', 'pooled_age_sex_models_top_human', 'pooled_age_sex_models_top_nh'), 
     file = "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/top-prelim-age-sex-pooled-RSFs.rda")
write.csv(pooled_age_sex_models_tab, 
          "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/prelim-age-sex-pooled-RSFs-AICoutputs.csv")


# ----------------------- FIT TOP AGE-SEX MODELS + COMBOS ---------------------- 
# competing (1) top human model, (2) top non-human/landscape model, and (3) a combination of the two
pooled_age_sex_models_final <- list()
pooled_age_sex_models_final[[1]] <- ksModel(age_sex, human6)
pooled_age_sex_models_final[[2]] <- ksModel(age_sex, nh11)
pooled_age_sex_models_final[[3]] <- ksModel(age_sex, pooled_combo1)
pooled_age_sex_models_final[[4]] <- ksModel(age_sex, pooled_combo2)
pooled_age_sex_models_final[[5]] <- ksModel(age_sex, pooled_combo3)
pooled_age_sex_models_final[[6]] <- ksModel(age_sex, null_mod)

names(pooled_age_sex_models_final) <- c('human6','nh11','pooled_combo1', 'pooled_combo2', 'pooled_combo3', 'null_mod')

# Save outputs
save(list = c("age_sex","pooled_age_sex_models_final"), 
     file = "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/final-age-sex-pooled-RSFs.rda")

# Create model selection table, exclude models w/ convergence issues
mod_select_final <- purrr::map(pooled_age_sex_models_final, function(x){!is.na(AIC(x)) & !is.nan(AIC(x))}) %>% unlist()
pooled_age_sex_models_final_good <- pooled_age_sex_models_final[mod_select_final] 
pooled_age_sex_models_final_tab <- bbmle::ICtab(pooled_age_sex_models_final_good, 
                                                mnames = names(pooled_age_sex_models_final_good), 
                                                type=c("AIC","BIC","AICc","qAIC","qAICc"),
                                                weights = FALSE, delta = TRUE, base = TRUE,
                                                logLik=FALSE, sort = TRUE,
                                                nobs=NULL, dispersion = 1, k = 2)
pooled_age_sex_models_final_tab <- as.data.frame(pooled_age_sex_models_final_tab)
# save top model
pooled_age_sex_models_final_top <- pooled_age_sex_models_final$pooled_combo3
summary(pooled_age_sex_models_final_top)
save(list = c('pooled_age_sex_models_final_tab', 'pooled_age_sex_models_final_top'), 
     file = "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/top-final-age-sex-pooled-RSF.rda")
write.csv(pooled_age_sex_models_final_tab, 
          "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/final-age-sex-pooled-RSFs-AICoutputs.csv")