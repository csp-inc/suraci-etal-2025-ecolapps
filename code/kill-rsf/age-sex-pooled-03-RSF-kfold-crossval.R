#install_github("BastilleRousseau/IndRSA")
library(devtools)
library(IndRSA)


# ------------------------------------------------------------------------------
# RUNNING K-FOLD CROSS-VALIDATION
# ------------------------------------------------------------------------------

# Load top model, perform cross validation
load(file = "model-outputs/top-final-age-sex-pooled-RSF.rda", verbose = TRUE) 

kfold_outputs_pooled <- kfoldRSF(
  pooled_age_sex_models_final_top,
  k = 5,
  nrepet = 50, # ran at 100 for a week without finishing
  nbins = 10,
  jitter = FALSE, # whether to add some random noise to the predictions (useful when the model is fitted on categorical variables, which can produces error in the ranking process)
  random = TRUE,
  method = method, # specified in function - spearman
  reproducible = TRUE # whether to use a fixed seed for each repetition
)
# export
write_csv(kfold_outputs_pooled, "cross-validation/kfold-cross-val-age-sex-pooled.csv")

# aggregate to get mean and SD
kfold_outputs_pooled_ag <- aggregate(. ~ type, kfold_outputs_pooled, function(x) c(mean = mean(x), sd = sd(x))) %>%
  mutate(mean = kfold[,"mean"],
         sd = kfold[,"sd"]) %>%
  dplyr::select(type, mean, sd)
# export
write_csv(kfold_outputs_pooled_ag, "cross-validation/mean-sd-kfold-cross-val-age-sex-pooled.csv")