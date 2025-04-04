## ---------------------------
##
## Script name: age-sex-03-RSF-kfold-crossval.R
##
## Purpose: Performs k-fold cross-validation on top RSF models for each age-sex
## class
##
## Author: Mae Lacey
##
## Email contact: mae@csp-inc.org
##
## ---------------------------


#install_github("BastilleRousseau/IndRSA")
library(devtools)
library(IndRSA)


# ------------------------------------------------------------------------------
# RUNNING K-FOLD CROSS-VALIDATION
# ------------------------------------------------------------------------------

# ADULT FEMALES: load top model, perform cross validation
load(file = "model-outputs/top-final-RSF-ad-fem.rda", verbose = TRUE) 

kfold_outputs_ad_fem <- kfoldRSF(
  age_sex_models_final_top,
  k = 5,
  nrepet = 100, # go back down to 50 if needed
  nbins = 10,
  jitter = FALSE, # whether to add some random noise to the predictions (useful when the model is fitted on categorical variables, which can produces error in the ranking process)
  random = TRUE,
  method = method, # specified in function - spearman
  reproducible = TRUE # whether to use a fixed seed for each repetition
)
# export
write_csv(kfold_outputs_ad_fem, "cross-validation/kfold-cross-val-ad-fem.csv")

# aggregate to get mean and SD
kfold_outputs_ad_fem_ag <- aggregate(. ~ type, kfold_outputs_ad_fem, function(x) c(mean = mean(x), sd = sd(x))) %>%
  mutate(mean = kfold[,"mean"],
         sd = kfold[,"sd"]) %>%
  dplyr::select(type, mean, sd)
# export
write_csv(kfold_outputs_ad_fem_ag, "cross-validation/mean-sd-kfold-cross-val-ad-fem.csv")


# ADULT MALES: load top model, perform cross validation
load(file = "model-outputs/top-final-RSF-ad-male.rda", verbose = TRUE) 

kfold_outputs_ad_male <- kfoldRSF(
  age_sex_models_final_top,
  k = 5,
  nrepet = 100, # go back down to 50 if needed
  nbins = 10,
  jitter = FALSE, 
  random = TRUE,
  method = method, 
  reproducible = TRUE 
)
# export
write_csv(kfold_outputs_ad_male, "cross-validation/kfold-cross-val-ad-male.csv")

# aggregate to get mean and SD
kfold_outputs_ad_male_ag <- aggregate(. ~ type, kfold_outputs_ad_male, function(x) c(mean = mean(x), sd = sd(x))) %>%
  mutate(mean = kfold[,"mean"],
         sd = kfold[,"sd"]) %>%
  dplyr::select(type, mean, sd)
# export
write_csv(kfold_outputs_ad_male_ag, "cross-validation/mean-sd-kfold-cross-val-ad-male.csv")


# DISPERSER FEMALES: load top model, perform cross validation
load(file = "model-outputs/top-final-RSF-disp-fem.rda", verbose = TRUE) 

kfold_outputs_disp_fem <- kfoldRSF(
  age_sex_models_final_top,
  k = 5,
  nrepet = 100, # go back down to 50 if needed
  nbins = 10,
  jitter = FALSE, 
  random = TRUE,
  method = method, 
  reproducible = TRUE 
)
# export
write_csv(kfold_outputs_disp_fem, "cross-validation/kfold-cross-val-disp-fem.csv")

# aggregate to get mean and SD
kfold_outputs_disp_fem_ag <- aggregate(. ~ type, kfold_outputs_disp_fem, function(x) c(mean = mean(x), sd = sd(x))) %>%
  mutate(mean = kfold[,"mean"],
         sd = kfold[,"sd"]) %>%
  dplyr::select(type, mean, sd)
# export
write_csv(kfold_outputs_disp_fem_ag, "cross-validation/mean-sd-kfold-cross-val-disp-fem.csv")


# DISPERSER MALES: load top model, perform cross validation
load(file = "model-outputs/top-final-RSF-disp-male.rda", verbose = TRUE) 

kfold_outputs_disp_male <- kfoldRSF(
  age_sex_models_final_top,
  k = 5,
  nrepet = 100, # go back down to 50 if needed
  nbins = 10,
  jitter = FALSE, 
  random = TRUE,
  method = method, 
  reproducible = TRUE 
)
# export
write_csv(kfold_outputs_disp_male, "cross-validation/kfold-cross-val-disp-male.csv")

# aggregate to get mean and SD
kfold_outputs_disp_male_ag <- aggregate(. ~ type, kfold_outputs_disp_male, function(x) c(mean = mean(x), sd = sd(x))) %>%
  mutate(mean = kfold[,"mean"],
         sd = kfold[,"sd"]) %>%
  dplyr::select(type, mean, sd)
# export
write_csv(kfold_outputs_disp_male_ag, "cross-validation/mean-sd-kfold-cross-val-disp-male.csv")