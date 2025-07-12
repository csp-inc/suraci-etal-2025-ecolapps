# Define models for age-sex class RSFs

# ------------------------------ MODEL FUNCTIONS -------------------------------
ksModel <- function(data, model){
  mod <- bglmer(model, data = data, family = binomial,
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                cov.prior = wishart, fixef.prior = NULL)
  return(mod)
}

ksModel_nm <- function(data, model){
  mod <- bglmer(model, data = data, family = binomial,
                control=glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=2e5)),
                cov.prior = wishart, fixef.prior = NULL)
  return(mod)
}


# ------------------------- MODEL SELECTION: STAGE 1 ---------------------------
# SUB-MODEL SELECTION

# human models

human1 <- formula(kill ~ agcov + devcov + roaddist + 
                    (1 + agcov + devcov + roaddist | Cougar_ID))

human2 <- formula(kill ~ agdist + roaddist + devdist + 
                    (1 + agdist + roaddist + devdist | Cougar_ID))

human3 <- formula(kill ~ agcov + devcov + roaddist + activity + 
                    (1 + agcov + devcov + roaddist + activity | Cougar_ID))

human4 <- formula(kill ~ agdist + roaddist + devdist + activity + 
                    (1 + agdist + roaddist + devdist + activity | Cougar_ID))

human5 <- formula(kill ~ agcov + agcov_sq + devcov + devcov_sq +  roaddist +  roaddist_sq +
                    (1 + agcov + devcov + roaddist | Cougar_ID))

human6 <- formula(kill ~ agdist + agdist_sq +  roaddist + roaddist_sq + devdist + devdist_sq + activity + 
                    (1 + agdist + roaddist + devdist + activity | Cougar_ID))

human7 <- formula(kill ~ agcov + agcov_sq + devcov + devcov_sq + 
                    (1 + agcov + devcov | Cougar_ID))

human8 <- formula(kill ~ agdist + agdist_sq +  devdist + devdist_sq + 
                    (1 + agdist + devdist | Cougar_ID))

human9 <- formula(kill ~ agcov + agcov_sq + devcov + devcov_sq + roaddist + 
                    (1 + agcov + devcov + roaddist | Cougar_ID))

human10 <- formula(kill ~ agdist + agdist_sq +  devdist + devdist_sq + activity + 
                     (1 + agdist + devdist + activity | Cougar_ID))

human11 <- formula(kill ~ agcov + agcov_sq + devcov + devcov_sq + activity + 
                     (1 + agcov + devcov + activity | Cougar_ID))

human12 <- formula(kill ~ agcov + agcov_sq + devcov + devcov_sq +  roaddist +  roaddist_sq + activity + 
                     (1 + agcov + devcov + roaddist + activity | Cougar_ID))

# non-human models

nh1 <- formula(kill ~ slope +  RAPtree + RAPtree_sq + shrubcov + shrubcov_sq + ripdist + 
                 (1 + slope + RAPtree + shrubcov +  ripdist | Cougar_ID)) # terrain_lcveg1

nh2 <- formula(kill ~ slope + slope_sq + edge + RAPtree + RAPtree_sq + shrubcov + shrubcov_sq + ripdist + 
                 (1 + slope + edge + RAPtree + shrubcov +  ripdist | Cougar_ID)) # terrain_lcveg2

nh3 <- formula(kill ~ slope + slope_sq + edge + RAPtree + RAPtree_sq + shrubcov + shrubcov_sq + 
                 (1 + slope + edge + RAPtree + shrubcov | Cougar_ID)) # terrain_lcveg3

nh4 <- formula(kill ~ slope + edge + RAPtree + shrubcov + 
                 (1 + slope + edge + RAPtree + shrubcov | Cougar_ID)) # terrain_lcveg4

nh5 <- formula(kill ~ edge + RAPtree + shrubcov + 
                 (1 + edge + RAPtree + shrubcov | Cougar_ID))

nh6 <- formula(kill ~ slope + edge + RAPtree + shrubcov + ripdist + 
                 (1 + slope + edge + RAPtree + shrubcov | Cougar_ID)) # terrain_lcveg_climate2

nh7 <- formula(kill ~ slope + slope_sq + edge + RAPtree + shrubcov +
                 (1 + slope + edge + RAPtree + shrubcov | Cougar_ID)) 

nh8 <- formula(kill ~ edge + RAPtree + RAPtree_sq + shrubcov + shrubcov_sq + 
                 (1 + edge + RAPtree + shrubcov | Cougar_ID))

nh9 <- formula(kill ~ RAPtree + RAPtree_sq + shrubcov + shrubcov_sq + ripdist +
                 (1 + RAPtree + shrubcov + ripdist | Cougar_ID))

nh10 <- formula(kill ~ edge + RAPtree + shrubcov + shrubcov_sq + ripdist + 
                  (1 + edge + RAPtree + shrubcov + ripdist | Cougar_ID))

nh11 <- formula(kill ~ slope + slope_sq + edge + RAPtree + RAPtree_sq + shrubcov + ripdist + 
                  (1 + slope + edge + RAPtree + shrubcov + ripdist | Cougar_ID)) # terrain_lcveg_climate3

nh12 <- formula(kill ~ slope + RAPtree + shrubcov + ripdist + 
                  (1 + slope + RAPtree + shrubcov + ripdist | Cougar_ID)) # terrain_lcveg_climate6

null_mod <- formula(kill ~ 1 + (1 | Cougar_ID))


# ------------------------- MODEL SELECTION: STAGE 2 ---------------------------

# INDIVIDUAL AGE/SEX CLASS MODELS: COMBINED MODEL SELECTION

# adult female
af_combo1 <- formula(kill ~ agdist + agdist_sq +  devdist + devdist_sq + slope + 
                       slope_sq + edge + RAPtree + RAPtree_sq + shrubcov + shrubcov_sq +
                       (1 + agdist + devdist + slope + edge + RAPtree + shrubcov | Cougar_ID))

af_combo2 <- formula(kill ~ agdist + agdist_sq +  devdist + devdist_sq + slope + 
                       slope_sq + edge + RAPtree + RAPtree_sq + shrubcov + shrubcov_sq + devdist*slope +
                       (1 + agdist + devdist + slope + edge + RAPtree + shrubcov | Cougar_ID))

# adult male
am_combo1 <- formula(kill ~ slope + slope_sq + edge + RAPtree + RAPtree_sq + shrubcov + 
                       agdist + agdist_sq +  devdist + devdist_sq + ripdist + 
                       (1 + slope + edge + RAPtree + shrubcov + ripdist + agdist + devdist | Cougar_ID)) 

am_combo2 <- formula(kill ~ slope + slope_sq + edge + RAPtree + RAPtree_sq + shrubcov + 
                       agdist + agdist_sq +  devdist + devdist_sq +
                       (1 + slope + edge + RAPtree + shrubcov + agdist + devdist | Cougar_ID)) 

am_combo3 <- formula(kill ~ slope + slope_sq + edge + RAPtree + RAPtree_sq + shrubcov + 
                       agdist + agdist_sq +  devdist + devdist_sq + devdist*slope +
                       (1 + slope + edge + RAPtree + shrubcov + agdist + devdist | Cougar_ID)) 

# disperser female
df_combo1 <- formula(kill ~ agdist + agdist_sq +  devdist + devdist_sq + 
                       RAPtree + RAPtree_sq + shrubcov + shrubcov_sq + ripdist +
                       (1 + RAPtree + shrubcov + ripdist + agdist + devdist | Cougar_ID))

df_combo2 <- formula(kill ~ agdist + agdist_sq +  devdist + devdist_sq + 
                       RAPtree + RAPtree_sq + shrubcov + shrubcov_sq +
                       (1 + RAPtree + shrubcov + agdist + devdist | Cougar_ID))

df_combo3 <- formula(kill ~ agdist + agdist_sq +  devdist + devdist_sq + 
                       RAPtree + RAPtree_sq + shrubcov + shrubcov_sq + slope + devdist*slope +
                       (1 + RAPtree + shrubcov + ripdist + agdist + devdist + slope | Cougar_ID))

# disperser male
dm_combo1 <- formula(kill ~ agdist + agdist_sq +  roaddist + roaddist_sq + devdist + devdist_sq + activity + 
                       slope + slope_sq + edge + RAPtree + RAPtree_sq + shrubcov + shrubcov_sq +
                       (1 + agdist + roaddist + devdist + activity + slope + edge + RAPtree + shrubcov| Cougar_ID))

dm_combo2 <- formula(kill ~ agdist + agdist_sq +  roaddist + roaddist_sq + devdist + devdist_sq + activity + 
                       slope + slope_sq + edge + RAPtree + RAPtree_sq + shrubcov + shrubcov_sq + devdist*slope +
                       (1 + agdist + roaddist + devdist + activity + slope + edge + RAPtree + shrubcov| Cougar_ID))

# POOLED AGE/SEX CLASS MODELS: COMBINED MODEL SELECTION
pooled_combo1 <- formula(kill ~ agdist + agdist_sq +  roaddist + roaddist_sq + devdist + devdist_sq + activity + 
                           slope + slope_sq + edge + RAPtree + RAPtree_sq + shrubcov + ripdist + 
                           (1 + agdist + roaddist + devdist + activity + slope + edge + RAPtree + shrubcov + ripdist | Cougar_ID))
  
pooled_combo2 <- formula(kill ~ agdist + agdist_sq +  roaddist + roaddist_sq + devdist + devdist_sq + activity + 
                           slope + slope_sq + edge + RAPtree + RAPtree_sq + shrubcov + 
                           (1 + agdist + roaddist + devdist + activity + slope + edge + RAPtree + shrubcov | Cougar_ID))

pooled_combo3 <- formula(kill ~ agdist + agdist_sq +  roaddist + roaddist_sq + devdist + devdist_sq + activity + 
                           slope + slope_sq + edge + RAPtree + RAPtree_sq + shrubcov + devdist*slope +
                           (1 + agdist + roaddist + devdist + activity + slope + edge + RAPtree + shrubcov | Cougar_ID))
