# Define models for SSF analyses

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# STAGE 1 - SUB-MODEL SELECTION

#------------- Define Models - HUMAN ------------------
human1 <- formula(case ~ -1 + ag + dev + d_road + log_sl + cos_ta + (1 | unique_step) + 
                    (0 + ag + dev + d_road | name))

human2 <- formula(case ~ -1 + d_ag + d_dev + d_road + log_sl + cos_ta + (1 | unique_step) + 
                    (0 + d_ag + d_dev + d_road | name))

human3 <- formula(case ~ -1 + ag + dev + d_road + forestry + log_sl + cos_ta + (1 | unique_step) + 
                    (0 + ag + dev + d_road + forestry | name))

human4 <- formula(case ~ -1 + d_ag + d_dev + d_road + forestry + log_sl + cos_ta + (1 | unique_step) + 
                    (0 + d_ag + d_dev + d_road + forestry | name))

human5 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + d_road + d_road_sq + log_sl + cos_ta + 
                    (1 | unique_step) + (0 + ag + dev + d_road | name))

human6 <- formula(case ~ -1 + d_ag + d_ag_sq + d_dev + d_dev_sq + d_road + d_road_sq + forestry + log_sl + cos_ta + 
                    (1 | unique_step) + (0 + d_ag + d_dev + d_road + forestry | name))

human7 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + log_sl + cos_ta + 
                    (1 | unique_step) + (0 + ag + dev | name))

human8 <- formula(case ~ -1 + d_ag + d_ag_sq + d_dev + d_dev_sq + log_sl + cos_ta + 
                    (1 | unique_step) + (0 + d_ag + d_dev | name))

human9 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + d_road + dev:night + log_sl + cos_ta + 
                    (1 | unique_step) + (0 + ag + dev + d_road | name))

human10 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + dev:night + log_sl + cos_ta + 
                     (1 | unique_step) + (0 + ag + dev | name))

#------------- Define Models - NON-HUMAN ------------------
nh1 <- formula(case ~ -1 + slope + tpi + tree + tree_sq + shrub + shrub_sq + d_riparian + log_sl + cos_ta + 
                 (1 | unique_step) + (0 + slope + tpi + tree + shrub + d_riparian | name))

nh2 <- formula(case ~ -1 + slope + slope_sq + edge + tree + tree_sq + shrub + shrub_sq + d_riparian + log_sl + cos_ta + 
                 (1 | unique_step) + (0 + slope + edge + tree + shrub + d_riparian | name))

nh3 <- formula(case ~ -1 + slope + slope_sq + edge + tree + tree_sq + shrub + shrub_sq + log_sl + cos_ta + 
                 (1 | unique_step) + (0 + slope + edge + tree + shrub | name))

nh4 <- formula(case ~ -1 + slope + edge + tree + shrub + log_sl + cos_ta + 
                 (1 | unique_step) + (0 + slope + edge + tree + shrub | name))

nh5 <- formula(case ~ -1 + tpi + edge + tree + shrub + log_sl + cos_ta + 
                 (1 | unique_step) + (0 + tpi + edge + tree + shrub | name))

nh6 <- formula(case ~ -1 + slope + tpi + edge + tree + shrub + d_riparian + log_sl + cos_ta + 
                 (1 | unique_step) + (0 + slope + tpi + edge + tree + shrub + d_riparian | name))

nh7 <- formula(case ~ -1 + slope + slope_sq + edge + tree + shrub + log_sl + cos_ta + 
                 (1 | unique_step) + (0 + slope + edge + tree + shrub | name))

nh8 <- formula(case ~ -1 + edge + tree + tree_sq + shrub + shrub_sq + log_sl + cos_ta + 
                 (1 | unique_step) + (0 + edge + tree + shrub | name))

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# STAGE 2 - COMBINED MODEL SELECTION

# Adult Female
af_combo1 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + d_road + d_road_sq + slope + slope_sq + 
                          edge + tree + tree_sq + shrub + shrub_sq + d_riparian + log_sl + cos_ta + 
                          (1 | unique_step) + (0 + ag + dev + d_road + slope + tree + shrub + d_riparian | name))
af_combo2 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + d_road + d_road_sq + slope + slope_sq + 
                          edge + tree + tree_sq + shrub + shrub_sq + log_sl + cos_ta + 
                          (1 | unique_step) + (0 + ag + dev + d_road + slope + tree + shrub | name))
af_combo3 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + d_road + d_road_sq + slope + slope_sq + 
                          edge + tree + tree_sq + shrub + shrub_sq + dev:slope + log_sl + cos_ta + 
                          (1 | unique_step) + (0 + ag + dev + d_road + slope + tree + shrub | name))


# Adult Male
am_combo1 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + d_road + d_road_sq + slope + slope_sq + 
                       edge + tree + tree_sq + shrub + shrub_sq + d_riparian + log_sl + cos_ta + 
                       (1 | unique_step) + (0 + ag + dev + d_road + slope + tree + shrub + d_riparian | name))
am_combo2 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + d_road + d_road_sq + slope + slope_sq + 
                       edge + tree + tree_sq + shrub + shrub_sq + d_riparian + dev:slope + log_sl + cos_ta + 
                       (1 | unique_step) + (0 + ag + dev + d_road + slope + tree + shrub + d_riparian | name))

# Disperser Female
df_combo1 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + slope + slope_sq + 
                       edge + tree + tree_sq + shrub + shrub_sq + d_riparian + log_sl + cos_ta + 
                       (1 | unique_step) + (0 + ag + dev + slope + tree + shrub + d_riparian | name))

df_combo2 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + slope + 
                       edge + tree + tree_sq + shrub + shrub_sq + d_riparian + log_sl + cos_ta + 
                       (1 | unique_step) + (0 + ag + dev + slope + tree + shrub + d_riparian | name))

df_combo3 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + slope + 
                       edge + tree + tree_sq + shrub + shrub_sq + d_riparian + dev:slope + log_sl + cos_ta + 
                       (1 | unique_step) + (0 + ag + dev + slope + tree + shrub + d_riparian | name))

# Disperser Male
dm_combo1 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + d_road + d_road_sq + slope + slope_sq + 
                       edge + tree + tree_sq + shrub + shrub_sq + d_riparian + log_sl + cos_ta + 
                       (1 | unique_step) + (0 + ag + dev + d_road + slope + tree + shrub + d_riparian | name))
dm_combo2 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + d_road + slope + slope_sq + 
                       edge + tree + tree_sq + d_riparian + log_sl + cos_ta + 
                       (1 | unique_step) + (0 + ag + dev + d_road + slope + tree + d_riparian | name))

# All Cats
ac_combo1 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + d_road + d_road_sq + slope + slope_sq + 
                       edge + tree + tree_sq + shrub + shrub_sq + d_riparian + log_sl + cos_ta + 
                       (1 | unique_step) + (0 + ag + dev + d_road + slope + tree + shrub + d_riparian | name))
ac_combo2 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + d_road + 
                       edge + tree + tree_sq + shrub + shrub_sq + log_sl + cos_ta + 
                       (1 | unique_step) + (0 + ag + dev + d_road + tree + shrub | name))
ac_combo3 <- formula(case ~ -1 + ag + ag_sq + dev + dev_sq + d_road + slope +
                       edge + tree + tree_sq + shrub + shrub_sq + dev:slope + log_sl + cos_ta + 
                       (1 | unique_step) + (0 + ag + dev + d_road + slope + tree + shrub | name))
