## ---------------------------
##
## Script name: 04-ssf-corss-val.R
##
## Purpose: Perform k-fold cross-validation
## on top SSF models for each age-sex class
##
## Author: Justin Suraci
##
## Email contact: justin@csp-inc.org
##
## ---------------------------

library(tidyverse)
library(glmmTMB)
source("utils/kfoldSSF_glmm.R")

# Load top models for each age/sex class
load("output/ad-fem-ssf-top-mod-v2.rda")
load("output/ad-male-ssf-top-mod-v2.rda")
load("output/disp-fem-ssf-top-mod-v2.rda")
load("output/disp-male-ssf-top-mod-v2.rda")
load("output/all-cats-ssf-top-mod.rda")

# Run k-fold cv for each model and save
ad_fem_cv <- kfoldSSF(mod = ad_fem_top, strata_name = "unique_step", k = 5, 
                      nrepet = 30, reproducible = FALSE, details = TRUE)

ad_male_cv <- kfoldSSF(mod = ad_male_top, strata_name = "unique_step", k = 5, 
                      nrepet = 30, reproducible = FALSE, details = TRUE)

disp_fem_cv <- kfoldSSF(mod = disp_fem_top, strata_name = "unique_step", k = 5, 
                      nrepet = 30, reproducible = FALSE, details = TRUE)

disp_male_cv <- kfoldSSF(mod = disp_male_top, strata_name = "unique_step", k = 5, 
                       nrepet = 30, reproducible = FALSE, details = TRUE)

all_cats_cv <- kfoldSSF(mod = all_cats_top, strata_name = "unique_step", k = 5, 
                         nrepet = 30, reproducible = FALSE, details = TRUE)

save(list = c("ad_fem_cv", "ad_male_cv", "disp_fem_cv", "disp_male_cv", "all_cats_cv"), file = "output/ssf-cross-val-30rep.rda")

cv_summarize <- function(dat, as_class){
  d <- dat %>% 
    filter(kfold != 0) %>% 
    group_by(type) %>% 
    summarise(mean = mean(kfold), sd = sd(kfold))
  
  d2 <- data.frame('age-sex class' = as_class, 
                   'Corr Observed' = paste0(round(d$mean[d$type == 'obs'], 2),
                                           " (", round(d$sd[d$type == 'obs'], 2), ")"),
                   'Corr Rand' = paste0(round(d$mean[d$type == 'rand'], 2),
                                        " (", round(d$sd[d$type == 'rand'], 2), ")"))
  return(d2)
}

out <- rbind(cv_summarize(ad_fem_cv, 'Adult female'),
      cv_summarize(ad_male_cv, 'Adult male'),
      cv_summarize(disp_fem_cv, 'Disperser female'),
      cv_summarize(disp_male_cv, 'Disperser male'),
      cv_summarize(all_cats_cv, 'All cats'))
colnames(out)<- c("Age/sex class", "Rho Observed", "Rho Random")
