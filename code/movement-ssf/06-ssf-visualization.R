## ---------------------------
##
## Script name: 06-ssf-visualization.R
##
## Purpose: Code for all SSF-related plots
##
## Author: Justin Suraci
##
## Email contact: justin@csp-inc.org
##
## ---------------------------

library(tidyverse)
library(terra)
library(sf)
library(glmmTMB)
# library(dotwhisker)
library(broom)
library(broom.mixed)
library(gridExtra)

# Get AOI
puma_crs <- st_crs(26910)
aoi = st_read("data/puma-refined-aoi.geojson") %>% st_transform(puma_crs)

# Get Analysis Datasets
load(file = "data/puma-ssf-dataset-20230701.rda", verbose = TRUE)
ssf_dat <- puma_ssf_dat %>% 
  dplyr::mutate(case = as.numeric(case_))

# Load top models for each age/sex class
load("output/ad-fem-ssf-top-mod-v2.rda")
load("output/ad-male-ssf-top-mod-v2.rda")
load("output/disp-fem-ssf-top-mod-v2.rda")
load("output/disp-male-ssf-top-mod-v2.rda")
load("output/all-cats-ssf-top-mod.rda")

# Set colors
af_col = '#7f32a8'
am_col = '#e3820b'
df_col = '#d184db'
dm_col = '#e3bf0b'
pooled_col = 'grey50'

#---------------------------------------
# COEF DOT AND WHISKER PLOTS
#---------------------------------------

# Prep model coefficients
af_coef <- broom.mixed::tidy(ad_fem_top) %>% filter(effect == 'fixed')
af_coef <- transform(af_coef,
                   term=sprintf("%s.%s", component, term))
am_coef <- broom.mixed::tidy(ad_male_top) %>% filter(effect == 'fixed')
am_coef <- transform(am_coef,
                   term=sprintf("%s.%s", component, term))
df_coef <- broom.mixed::tidy(disp_fem_top) %>% filter(effect == 'fixed')
df_coef <- transform(df_coef,
                     term=sprintf("%s.%s", component, term))
dm_coef <- broom.mixed::tidy(disp_male_top) %>% filter(effect == 'fixed')
dm_coef <- transform(dm_coef,
                     term=sprintf("%s.%s", component, term))
all_coef <- broom.mixed::tidy(all_cats_top) %>% filter(effect == 'fixed')
all_coef <- transform(all_coef,
                     term=sprintf("%s.%s", component, term))


# Adult Female
pdf(file = "output/adult-female-ssf-coefs.pdf", width = 5, height = 4)
af_coef_plot <- dwplot(af_coef,
                     ci = 0.95,
                     dot_args = list(size = 2.5, shape = 21, fill = af_col, color = 'black'),
                     whisker_args = list(size = 1, color = "black"),
                     vline = geom_vline(xintercept = 0,
                                        colour = "grey60",
                                        linetype = 2)) %>% 
  relabel_predictors(cond.edge = "Annual edge % cover", 
                     cond.tree = "Annual tree % cover",
                     cond.tree_sq = "Annual tree^2",
                     cond.shrub = "Shrub % cover",
                     cond.shrub_sq = "Shrub^2",
                     cond.ag = "Ag % cover",
                     cond.ag_sq = "Ag^2",
                     cond.dev = "Developed % cover",
                     cond.dev_sq = "Developed^2",
                     cond.d_road = "Dist to road",
                     cond.d_road_sq = "Dist to road^2",
                     cond.slope = "Slope",
                     cond.slope_sq = "Slope^2",
                     cond.log_sl = "[Log step length]",
                     cond.cos_ta = "[Cosine turn angle]"
  ) +
  theme_bw() + 
  xlab("") + ylab("") + 
  theme(legend.position = "none")
af_coef_plot
dev.off()

# Adult Male
pdf(file = "output/adult-male-ssf-coefs.pdf", width = 5, height = 4)
am_coef_plot <- dwplot(am_coef,
                       ci = 0.95,
                       dot_args = list(size = 2.5, shape = 21, fill = am_col, color = 'black'),
                       whisker_args = list(size = 1, color = "black"),
                       vline = geom_vline(xintercept = 0,
                                          colour = "grey60",
                                          linetype = 2)) %>% 
  relabel_predictors(cond.edge = "Annual edge % cover", 
                     cond.tree = "Annual tree % cover",
                     cond.tree_sq = "Annual tree^2",
                     cond.shrub = "Shrub % cover",
                     cond.shrub_sq = "Shrub^2",
                     cond.ag = "Ag % cover",
                     cond.ag_sq = "Ag^2",
                     cond.dev = "Developed % cover",
                     cond.dev_sq = "Developed^2",
                     cond.d_road = "Dist to road",
                     cond.d_road_sq = "Dist to road^2",
                     cond.slope = "Slope",
                     cond.slope_sq = "Slope^2",
                     cond.d_riparian = "Dist to riparian",
                     'cond.dev:slope' = "Developed x Slope",
                     cond.log_sl = "[Log step length]",
                     cond.cos_ta = "[Cosine turn angle]"
  ) +
  theme_bw() + 
  xlab("") + ylab("") + 
  theme(legend.position = "none")
am_coef_plot
dev.off()

# Disperser Female
pdf(file = "output/disperser-female-ssf-coefs.pdf", width = 5, height = 4)
df_coef_plot <- dwplot(df_coef,
                       ci = 0.95,
                       dot_args = list(size = 2.5, shape = 21, fill = df_col, color = 'black'),
                       whisker_args = list(size = 1, color = "black"),
                       vline = geom_vline(xintercept = 0,
                                          colour = "grey60",
                                          linetype = 2)) %>% 
  relabel_predictors(cond.edge = "Annual edge % cover", 
                     cond.tree = "Annual tree % cover",
                     cond.tree_sq = "Annual tree^2",
                     cond.shrub = "Shrub % cover",
                     cond.shrub_sq = "Shrub^2",
                     cond.ag = "Ag % cover",
                     cond.ag_sq = "Ag^2",
                     cond.dev = "Developed % cover",
                     cond.dev_sq = "Developed^2",
                     cond.slope = "Slope",
                     cond.d_riparian = "Dist to riparian",
                     cond.log_sl = "[Log step length]",
                     cond.cos_ta = "[Cosine turn angle]"
  ) +
  theme_bw() + 
  xlab("") + ylab("") + 
  theme(legend.position = "none")
df_coef_plot
dev.off()

# Disperser Male
pdf(file = "output/disperser-male-ssf-coefs.pdf", width = 5, height = 4)
dm_coef_plot <- dwplot(dm_coef,
                       ci = 0.95,
                       dot_args = list(size = 2.5, shape = 21, fill = dm_col, color = 'black'),
                       whisker_args = list(size = 1, color = "black"),
                       vline = geom_vline(xintercept = 0,
                                          colour = "grey60",
                                          linetype = 2)) %>% 
  relabel_predictors(cond.edge = "Annual edge % cover", 
                     cond.tree = "Annual tree % cover",
                     cond.tree_sq = "Annual tree^2",
                     cond.ag = "Ag % cover",
                     cond.ag_sq = "Ag^2",
                     cond.dev = "Developed % cover",
                     cond.dev_sq = "Developed^2",
                     cond.d_road = "Dist to road",
                     cond.slope = "Slope",
                     cond.slope_sq = "Slope^2",
                     cond.d_riparian = "Dist to riparian",
                     cond.log_sl = "[Log step length]",
                     cond.cos_ta = "[Cosine turn angle]"
  ) +
  theme_bw() + 
  xlab("") + ylab("") + 
  theme(legend.position = "none")
dm_coef_plot
dev.off()

# ALL AGE/SEX CLASSES IN FOUR PANELS
pdf(file = "output/all-age-sex-ssf-coefs.pdf", width = 15, height = 4)
grid.arrange(af_coef_plot, am_coef_plot, df_coef_plot, dm_coef_plot, nrow = 1)
dev.off()
 
#------------
# ALL CATS POOLED
pdf(file = "output/pooled-ssf-coefs.pdf", width = 5, height = 4)
all_coef_plot <- dwplot(all_coef,
                       ci = 0.95,
                       dot_args = list(size = 2.5, shape = 21, fill = pooled_col, color = 'black'),
                       whisker_args = list(size = 1, color = "black"),
                       vline = geom_vline(xintercept = 0,
                                          colour = "grey60",
                                          linetype = 2)) %>% 
  relabel_predictors(cond.edge = "Annual edge % cover", 
                     cond.tree = "Annual tree % cover",
                     cond.tree_sq = "Annual tree^2",
                     cond.shrub = "Shrub % cover",
                     cond.shrub_sq = "Shrub^2",
                     cond.ag = "Ag % cover",
                     cond.ag_sq = "Ag^2",
                     cond.dev = "Developed % cover",
                     cond.dev_sq = "Developed^2",
                     cond.d_road = "Dist to road",
                     cond.d_road_sq = "Dist to road^2",
                     cond.slope = "Slope",
                     cond.slope_sq = "Slope^2",
                     cond.d_riparian = "Dist to riparian",
                     cond.log_sl = "[Log step length]",
                     cond.cos_ta = "[Cosine turn angle]"
  ) +
  theme_bw() + 
  xlab("") + ylab("") + 
  theme(legend.position = "none")
all_coef_plot
dev.off()


#---------------------------------------
# PLOT NON-LINEAR COVARIATE EFFECTS
#---------------------------------------

# FUNCTIONS

# Prep age/sex class specific datasets
# (recycled from analysis scripts)
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

# Linear prediction function
modPredict <- function(covs, mod, scale = "exp", keepFirst = FALSE, constrain = NULL){
  mod <- summary(mod)
  coefs <- mod$coefficients$cond[,1]
  se <- mod$coefficients$cond[,2]
  if(keepFirst == FALSE){
    coefs <- coefs[-1]
    se <- se[-1]
  }
  
  # Create linear prediction (and transform, as specified in function call)
  pred <- (covs %*% coefs)
  se_up <- (covs %*% (coefs + se))
  se_low <- (covs %*% (coefs - se))
  if(!is.null(constrain)) pred <- pConstrain(pred, constrain)
  if(scale == "exp"){
    pred <- exp(pred)
    se_up <- exp(se_up)
    se_low <- exp(se_low)
  } 
  
  if(scale == "response"){
    pred <- exp(pred)/(1+exp(pred))
    se_up <- exp(se_up)/(1+exp(se_up))
    se_low <- exp(se_low)/(1+exp(se_low))
  } 
  pred_out <- data.frame(pred = pred, se_up = se_up,  se_low = se_low)
  return(pred_out)
}

# Unscale covs
unscale <- function(new_cov, orig_name, data){
  fcov = data %>% pull(orig_name)
  m_cov = mean(fcov, na.rm = T)
  sd_cov = sd(fcov, na.rm = T)
  
  out = (new_cov * sd_cov) + m_cov
  return(out)
}

# Create non-linear plot for selected covariate (and its polynomial term)
# mod: model to predict from
# cov_name: name (in quotes) of covariate. Function requires that the column name of the
# squared term for the selected cov takes the from "{cov_name}_sq"
# orig_name: original name of the covariate prior to data processing 
# data: data set used to fit the model specified in the mod argument
# x_name: x axis label on plot
# y_name: y axis label on plot
nonLinPlot <- function(mod, cov_name, orig_name, data, x_name = "", y_name = "", x_round = 0, x_div = NA, 
                       se_col = "#92d7e7", m_col = 'grey60', scale = "exp", keepFirst = FALSE, constrain = NULL){
  # Prep empty dataframe
  if(keepFirst) d_names <- names(fixef(mod)$cond) else d_names <- names(fixef(mod)$cond[-1])
  d_new <- array(0, dim = c(100, length(d_names)))
  colnames(d_new) <- d_names
  
  # Fill empty dataframe
  minMax <- data %>% pull(cov_name) %>% quantile(probs = c(0,1), na.rm = T)
  cov_new <- seq(minMax[1], minMax[2], length.out = 100)
  cov_sq_new <- cov_new^2
  d_cov_new <- d_new
  d_cov_new[,which(d_names == cov_name)] <- cov_new
  d_cov_new[,which(d_names == paste0(cov_name, "_sq"))] <- cov_sq_new
  d_cov_pred <- modPredict(as.matrix(d_cov_new), mod = mod, scale = scale, keepFirst = keepFirst, constrain = constrain)
  
  if(is.na(x_div)){
    outPlot <- ggplot(data = d_cov_pred, aes(x = cov_new, y = pred, ymin = se_low, ymax = se_up)) + 
      geom_ribbon(alpha = 0.4, fill = se_col) + 
      geom_line(color = m_col) +
      theme_bw() +
      xlab(x_name) + 
      ylab(y_name) +
      scale_x_continuous(breaks = seq(minMax[1], minMax[2], length.out = 5), 
                         labels = unscale(seq(minMax[1], minMax[2], length.out = 5), orig_name, data) %>% round(x_round))
  }
  
  if(!is.na(x_div)){
    outPlot <- ggplot(data = d_cov_pred, aes(x = cov_new, y = pred, ymin = se_low, ymax = se_up)) + 
      geom_ribbon(alpha = 0.4, fill = se_col) + 
      geom_line(color = m_col) +
      theme_bw() +
      xlab(x_name) + 
      ylab(y_name) +
      scale_x_continuous(breaks = seq(minMax[1], minMax[2], length.out = 5), 
                         labels = round(unscale(seq(minMax[1], minMax[2], length.out = 5), orig_name, data)/x_div, x_round))
  }
  
  out<-list(outPlot,d_cov_pred)
  names(out)<-c('plot','data')
  return(out)
}

#---------------------
# RUN FOR POOLED MOD
# prep dataset
pooled <- ssf_dat %>% datPrep()

# tree
p_tree <- nonLinPlot(mod = all_cats_top,
                      cov_name = "tree",
                      orig_name = "trCover_pcov_1km",
                      data = pooled,
                      y_name = "Relative probability of selection",
                      x_name = "Annual tree % cover",
                      keepFirst = TRUE) 
pdf(file = "output/pooled-tree-nonlin.pdf", width = 5, height = 5)
p_tree$plot
dev.off()

p_shrub <- nonLinPlot(mod = all_cats_top,
                     cov_name = "shrub",
                     orig_name = "shrub_pcov_150m",
                     data = pooled,
                     x_div = 0.01,
                     y_name = "Relative probability of selection",
                     x_name = "Shrub % cover",
                     keepFirst = TRUE)
pdf(file = "output/pooled-shrub-nonlin.pdf", width = 5, height = 5)
p_shrub$plot
dev.off()

p_ag <- nonLinPlot(mod = all_cats_top,
                      cov_name = "ag",
                      orig_name = "ag_pcov_150m",
                      data = pooled,
                      x_div = 0.01,
                      # x_round = 2,
                      x_name = "Agriculture % cover",
                      y_name = "Relative probability of selection",
                      keepFirst = TRUE)
pdf(file = "output/pooled-ag-nonlin.pdf", width = 5, height = 5)
p_ag$plot
dev.off()

p_dev <- nonLinPlot(mod = all_cats_top,
                      cov_name = "dev",
                      orig_name = "developed_150m",
                      data = pooled,
                      x_div = 0.01,
                      # x_round = 2,
                      x_name = "Development % cover",
                      y_name = "Relative probability of selection",
                      keepFirst = TRUE)
pdf(file = "output/pooled-dev-nonlin.pdf", width = 5, height = 5)
p_dev$plot
dev.off()

p_road <- nonLinPlot(mod = all_cats_top,
                    cov_name = "d_road",
                    orig_name = "road_dist",
                    data = pooled,
                    x_div = 1000,
                    # x_round = 2,
                    x_name = "Distance to road (km)",
                    y_name = "Relative probability of selection",
                    keepFirst = TRUE)
pdf(file = "output/pooled-road-nonlin.pdf", width = 5, height = 5)
p_road$plot
dev.off()

p_slope <- nonLinPlot(mod = all_cats_top,
                     cov_name = "slope",
                     orig_name = "slope_150m_path",
                     data = pooled,
                     # x_round = 2,
                     x_name = "Slope (degrees)",
                     y_name = "Relative probability of selection",
                     keepFirst = TRUE)
pdf(file = "output/pooled-slope-nonlin.pdf", width = 5, height = 5)
p_slope$plot
dev.off()


#---------------------
#---------------------
#RUN FOR AGE SEX CLASSES
# prep datasets
ad_fem <- ssf_dat %>% filter(age_class == "adult", sex == "Female") %>% datPrep()
ad_male <- ssf_dat %>% filter(age_class == "adult", sex == "Male") %>% datPrep()
disp_fem <- ssf_dat %>% filter(age_class == "disperser", sex == "Female") %>% datPrep()
disp_male <- ssf_dat %>% filter(age_class == "disperser", sex == "Male") %>% datPrep()

# tree
af_tree <- nonLinPlot(mod = ad_fem_top,
                      cov_name = "tree",
                      orig_name = "trCover_pcov_1km",
                      data = ad_fem,
                      y_name = "Relative probability of selection",
                      keepFirst = TRUE,
                      m_col = af_col,
                      se_col = af_col)

am_tree <- nonLinPlot(mod = ad_male_top,
                      cov_name = "tree",
                      orig_name = "trCover_pcov_1km",
                      data = ad_male,
                      keepFirst = TRUE,
                      m_col = am_col,
                      se_col = am_col)

df_tree <- nonLinPlot(mod = disp_fem_top,
                      cov_name = "tree",
                      orig_name = "trCover_pcov_1km",
                      data = disp_fem,
                      x_name = "Annual Tree % cover",
                      y_name = "Relative probability of selection",
                      keepFirst = TRUE,
                      m_col = df_col,
                      se_col = df_col)

dm_tree <- nonLinPlot(mod = disp_male_top,
                      cov_name = "tree",
                      orig_name = "trCover_pcov_1km",
                      data = disp_male,
                      x_name = "Annual Tree % cover",
                      keepFirst = TRUE,
                      m_col = dm_col,
                      se_col = dm_col)
pdf(file = 'output/tree-nonlin-plot.pdf', height = 6, width = 6)
grid.arrange(af_tree$plot, am_tree$plot, df_tree$plot, dm_tree$plot, nrow = 2)
dev.off()


# dev
af_dev <- nonLinPlot(mod = ad_fem_top,
                      cov_name = "dev",
                      orig_name = "developed_150m",
                      data = ad_fem,
                      y_name = "Relative probability of selection",
                      keepFirst = TRUE,
                      x_round = 2,
                      m_col = af_col,
                      se_col = af_col)

am_dev <- nonLinPlot(mod = ad_male_top,
                      cov_name = "dev",
                      orig_name = "developed_150m",
                      data = ad_male,
                      keepFirst = TRUE,
                      x_round = 2,
                      m_col = am_col,
                      se_col = am_col)

df_dev <- nonLinPlot(mod = disp_fem_top,
                      cov_name = "dev",
                      orig_name = "developed_150m",
                      data = disp_fem,
                      x_name = "Developed % cover",
                      y_name = "Relative probability of selection",
                      keepFirst = TRUE,
                      x_round = 2,
                      m_col = df_col,
                      se_col = df_col)

dm_dev <- nonLinPlot(mod = disp_male_top,
                      cov_name = "dev",
                      orig_name = "developed_150m",
                      data = disp_male,
                      x_name = "Developed % cover",
                      keepFirst = TRUE,
                      x_round = 2,
                      m_col = dm_col,
                      se_col = dm_col)
pdf(file = 'output/dev-nonlin-plot.pdf', height = 6, width = 6)
grid.arrange(af_dev$plot, am_dev$plot, df_dev$plot, dm_dev$plot, nrow = 2)
dev.off()

# ag
af_ag <- nonLinPlot(mod = ad_fem_top,
                     cov_name = "ag",
                     orig_name = "ag_pcov_150m",
                     data = ad_fem,
                     y_name = "Relative probability of selection",
                     keepFirst = TRUE,
                     x_round = 2,
                     m_col = af_col,
                     se_col = af_col)

am_ag <- nonLinPlot(mod = ad_male_top,
                     cov_name = "ag",
                     orig_name = "ag_pcov_150m",
                     data = ad_male,
                     keepFirst = TRUE,
                     x_round = 2,
                     m_col = am_col,
                     se_col = am_col)

df_ag <- nonLinPlot(mod = disp_fem_top,
                     cov_name = "ag",
                     orig_name = "ag_pcov_150m",
                     data = disp_fem,
                     x_name = "Agriculture % cover",
                     y_name = "Relative probability of selection",
                     keepFirst = TRUE,
                     x_round = 2,
                     m_col = df_col,
                     se_col = df_col)

dm_ag <- nonLinPlot(mod = disp_male_top,
                     cov_name = "ag",
                     orig_name = "ag_pcov_150m",
                     data = disp_male,
                     x_name = "Agriculture % cover",
                     keepFirst = TRUE,
                     x_round = 2,
                     m_col = dm_col,
                     se_col = dm_col)
pdf(file = 'output/ag-nonlin-plot.pdf', height = 6, width = 6)
grid.arrange(af_ag$plot, am_ag$plot, df_ag$plot, dm_ag$plot, nrow = 2)
dev.off()

# Road
am_road <- nonLinPlot(mod = ad_male_top,
                    cov_name = "d_road",
                    orig_name = "road_dist",
                    data = ad_male,
                    x_name = "Distance to road",
                    keepFirst = TRUE,
                    x_round = -2,
                    m_col = am_col,
                    se_col = am_col)

af_road <- nonLinPlot(mod = ad_fem_top,
                      cov_name = "d_road",
                      orig_name = "road_dist",
                      data = ad_fem,
                      y_name = "Relative probability of selection",
                      x_name = "Distance to road",
                      keepFirst = TRUE,
                      x_round = -2,
                      m_col = af_col,
                      se_col = af_col)
pdf(file = 'output/road-nonlin-plot-af-am.pdf', height = 3, width = 6)
grid.arrange(af_road$plot, am_road$plot, nrow = 1)
dev.off()

# slope
af_slope <- nonLinPlot(mod = ad_fem_top,
                    cov_name = "slope",
                    orig_name = "slope_150m_path",
                    data = ad_fem,
                    y_name = "Relative probability of selection",
                    keepFirst = TRUE,
                    x_round = 0,
                    m_col = af_col,
                    se_col = af_col)

am_slope <- nonLinPlot(mod = ad_male_top,
                    cov_name = "slope",
                    orig_name = "slope_150m_path",
                    x_name = "Slope",
                    data = ad_male,
                    keepFirst = TRUE,
                    x_round = 0,
                    m_col = am_col,
                    se_col = am_col)

dm_slope <- nonLinPlot(mod = disp_male_top,
                    cov_name = "slope",
                    orig_name = "slope_150m_path",
                    data = disp_male,
                    keepFirst = TRUE,
                    x_round = 0,
                    m_col = dm_col,
                    se_col = dm_col)

grid.arrange(af_slope$plot, am_slope$plot, dm_slope$plot, nrow = 1)

# Shrub
af_shrub <- nonLinPlot(mod = ad_fem_top,
                       cov_name = "shrub",
                       orig_name = "shrub_pcov_150m",
                       data = ad_fem,
                       y_name = "Relative probability of selection",
                       keepFirst = TRUE,
                       x_round = 2,
                       m_col = af_col,
                       se_col = af_col)

am_shrub <- nonLinPlot(mod = ad_male_top,
                       cov_name = "shrub",
                       orig_name = "shrub_pcov_150m",
                       x_name = "shrub",
                       data = ad_male,
                       keepFirst = TRUE,
                       x_round = 2,
                       m_col = am_col,
                       se_col = am_col)

df_shrub <- nonLinPlot(mod = disp_fem_top,
                       cov_name = "shrub",
                       orig_name = "shrub_pcov_150m",
                       data = disp_fem,
                       keepFirst = TRUE,
                       x_round = 2,
                       m_col = df_col,
                       se_col = df_col)

grid.arrange(af_shrub$plot, am_shrub$plot, df_shrub$plot, nrow = 1)
