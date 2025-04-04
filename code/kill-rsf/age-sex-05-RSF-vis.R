## ---------------------------
##
## Script name: age-sex-05-RSF-vis.R
##
## Purpose: Code for all RSF-related plots
##
## Author: Mae Lacey
##
## Email contact: mae@csp-inc.org
##
## ---------------------------


# ------------------------------------------------------------------------------ 
# ------------------------------ Coefficient plots -----------------------------
# ------------------------------------------------------------------------------ 

library(dotwhisker)
library(broom)
library(broom.mixed)
library(gridExtra)
library(dplyr)

# Load final models
load(file = "age-sex-class-RSFs/model-outputs/top-final-RSF-ad-fem.rda", verbose = TRUE) # load RSFs
ad_fem_models_final_top = age_sex_models_final_top
load(file = "age-sex-class-RSFs/model-outputs/top-final-RSF-ad-male.rda", verbose = TRUE) # load RSFs
ad_male_models_final_top = age_sex_models_final_top
load(file = "age-sex-class-RSFs/model-outputs/top-final-RSF-disp-fem.rda", verbose = TRUE) # load RSFs
disp_fem_models_final_top = age_sex_models_final_top
load(file = "age-sex-class-RSFs/model-outputs/top-final-RSF-disp-male.rda", verbose = TRUE) # load RSFs
disp_male_models_final_top = age_sex_models_final_top


# ------------------------------------------------------------------------------ 
# --------------------------- Model coefficient plots --------------------------
# ------------------------------------------------------------------------------ 

# Prep model coefficients
adfem_coef <- broom.mixed::tidy(ad_fem_models_final_top) %>% filter(effect == 'fixed')
adfem_coef <- transform(adfem_coef,
                        term=sprintf("%s.%s", effect, term))[-1,]
admale_coef <- broom.mixed::tidy(ad_male_models_final_top) %>% filter(effect == 'fixed')
admale_coef <- transform(admale_coef,
                         term=sprintf("%s.%s", effect, term))[-1,]
dispfem_coef <- broom.mixed::tidy(disp_fem_models_final_top) %>% filter(effect == 'fixed')
dispfem_coef <- transform(dispfem_coef,
                          term=sprintf("%s.%s", effect, term))[-1,]
dispmale_coef <- broom.mixed::tidy(disp_male_models_final_top) %>% filter(effect == 'fixed')
dispmale_coef <- transform(dispmale_coef,
                           term=sprintf("%s.%s", effect, term))[-1,]

# save CSVs for later
write_csv(adfem_coef, "age-sex-class-RSFs/model-outputs/adfem_coef.csv")
write_csv(admale_coef, "age-sex-class-RSFs/model-outputs/admale_coef.csv")
write_csv(dispfem_coef, "age-sex-class-RSFs/model-outputs/dispfem_coef.csv")
write_csv(dispmale_coef, "age-sex-class-RSFs/model-outputs/dispmale_coef.csv")


# ------------------------------ Individual plots ------------------------------

agesex <- c("#7f32a8", "#e3820b", "#d184db", "#e3bf0b")

# Adult Female
pdf(file = "age-sex-class-RSFs/model-outputs/age-sex-adult-female-rsf-coefs.pdf", width = 5, height = 4)
adfem_coef_plot <- dwplot(adfem_coef,
                          ci = 0.95,
                          dot_args = list(size = 2.5, shape = 22, fill = agesex[1], color = 'black'),
                          whisker_args = list(size = 1, color = "black"),
                          vline = geom_vline(xintercept = 0,
                                             colour = "grey60",
                                             linetype = 2)) %>% 
  relabel_predictors(fixed.edge = "Annual edge % cover",
                     fixed.RAPtree = "Annual tree % cover",
                     fixed.RAPtree_sq = "Annual tree^2",
                     fixed.shrubcov = "Shrub % cover",
                     fixed.shrubcov_sq = "Shrub^2",
                     fixed.agdist = "Dist to ag",
                     fixed.agdist_sq = "Dist to ag^2",
                     fixed.devdist = "Dist to developed",
                     fixed.devdist_sq = "Dist to developed^2",
                     fixed.roaddist = "Dist to road",
                     fixed.roaddist_sq = "Dist to road^2",
                     fixed.slope = "Slope",
                     fixed.slope_sq = "Slope^2",
                     fixed.activity = "Forestry activity",
                     fixed.ripdist = "Dist to riparian",
                     `fixed.devdist:slope` = "Dist to developed * Slope"
  ) + 
  theme_bw() + 
  xlab("") + ylab("") + 
  theme(legend.position = "none") 
adfem_coef_plot
dev.off()


# Adult Male
pdf(file = "age-sex-class-RSFs/model-outputs/age-sex-adult-male-rsf-coefs.pdf", width = 5, height = 4)
admale_coef_plot <- dwplot(admale_coef,
                           ci = 0.95,
                           dot_args = list(size = 2.5, shape = 22, fill = agesex[2], color = 'black'),
                           whisker_args = list(size = 1, color = "black"),
                           vline = geom_vline(xintercept = 0,
                                              colour = "grey60",
                                              linetype = 2)) %>% 
  relabel_predictors(fixed.edge = "Annual edge % cover",
                     fixed.RAPtree = "Annual tree % cover",
                     fixed.RAPtree_sq = "Annual tree^2",
                     fixed.shrubcov = "Shrub % cover",
                     fixed.shrubcov_sq = "Shrub^2",
                     fixed.agdist = "Dist to ag",
                     fixed.agdist_sq = "Dist to ag^2",
                     fixed.devdist = "Dist to developed",
                     fixed.devdist_sq = "Dist to developed^2",
                     fixed.roaddist = "Dist to road",
                     fixed.roaddist_sq = "Dist to road^2",
                     fixed.slope = "Slope",
                     fixed.slope_sq = "Slope^2",
                     fixed.activity = "Forestry activity",
                     fixed.ripdist = "Dist to riparian",
                     `fixed.devdist:slope` = "Dist to developed * Slope"
  ) + 
  theme_bw() + 
  xlab("") + ylab("") + 
  theme(legend.position = "none") 
admale_coef_plot
dev.off()

# Disperser Female
pdf(file = "age-sex-class-RSFs/model-outputs/age-sex-disperser-female-rsf-coefs.pdf", width = 5, height = 4)
dispfem_coef_plot <- dwplot(dispfem_coef,
                            ci = 0.95,
                            dot_args = list(size = 2.5, shape = 22, fill = agesex[3], color = 'black'),
                            whisker_args = list(size = 1, color = "black"),
                            vline = geom_vline(xintercept = 0,
                                               colour = "grey60",
                                               linetype = 2)) %>% 
  relabel_predictors(fixed.edge = "Annual edge % cover",
                     fixed.RAPtree = "Annual tree % cover",
                     fixed.RAPtree_sq = "Annual tree^2",
                     fixed.shrubcov = "Shrub % cover",
                     fixed.shrubcov_sq = "Shrub^2",
                     fixed.agdist = "Dist to ag",
                     fixed.agdist_sq = "Dist to ag^2",
                     fixed.devdist = "Dist to developed",
                     fixed.devdist_sq = "Dist to developed^2",
                     fixed.roaddist = "Dist to road",
                     fixed.roaddist_sq = "Dist to road^2",
                     fixed.slope = "Slope",
                     fixed.slope_sq = "Slope^2",
                     fixed.activity = "Forestry activity",
                     fixed.ripdist = "Dist to riparian",
                     `fixed.devdist:slope` = "Dist to developed * Slope"
  ) + 
  theme_bw() + 
  xlab("") + ylab("") + 
  theme(legend.position = "none") 
dispfem_coef_plot
dev.off()

# Disperser Male
pdf(file = "age-sex-class-RSFs/model-outputs/age-sex-disperser-male-rsf-coefs.pdf", width = 5, height = 4)
dispmale_coef_plot <- dwplot(dispmale_coef,
                             ci = 0.95,
                             dot_args = list(size = 2.5, shape = 22, fill = agesex[4], color = 'black'),
                             whisker_args = list(size = 1, color = "black"),
                             vline = geom_vline(xintercept = 0,
                                                colour = "grey60",
                                                linetype = 2)) %>% 
  relabel_predictors(fixed.edge = "Annual edge % cover",
                     fixed.RAPtree = "Annual tree % cover",
                     fixed.RAPtree_sq = "Annual tree^2",
                     fixed.shrubcov = "Shrub % cover",
                     fixed.shrubcov_sq = "Shrub^2",
                     fixed.agdist = "Dist to ag",
                     fixed.agdist_sq = "Dist to ag^2",
                     fixed.devdist = "Dist to developed",
                     fixed.devdist_sq = "Dist to developed^2",
                     fixed.roaddist = "Dist to roads",
                     fixed.roaddist_sq = "Dist to roads^2",
                     fixed.slope = "Slope",
                     fixed.slope_sq = "Slope^2",
                     fixed.activity = "Forestry activity",
                     fixed.ripdist = "Dist to riparian",
                     `fixed.devdist:slope` = "Dist to developed * Slope"
  ) + 
  theme_bw() + 
  xlab("") + ylab("") + 
  theme(legend.position = "none") 
dispmale_coef_plot
dev.off()


# ----------------------------- All plots combined -----------------------------
pdf(file = "age-sex-class-RSFs/model-outputs/all-age-sex-rsf-coefs.pdf", width = 15, height = 4)
grid.arrange(adfem_coef_plot, admale_coef_plot, dispfem_coef_plot, dispmale_coef_plot, nrow = 1)
dev.off()


# ------------------------------------------------------------------------------ 
# ---------------------- Plot non-linear covariate effects ---------------------
# ------------------------------------------------------------------------------ 

# Load raw data table
load("age-sex-class-RSFs/age_sex_RSF_covariate_extract_all.rda")

# --------------------------------- Functions ---------------------------------

# Prep age/sex class specific datasets
datPrep <- function(dat){
  dout <- dat %>% dplyr::select(Cougar_ID, Report_Id, age_sex, year, kill, ag_dist,
                                all_ag_pcov_100m, aspect_100m, developed_dist, 
                                developed_pcov_1km, FPA_TH_activity, riparian_dist, 
                                road_dist, shrub_pcov_500m, slope_1km, RAP_treecover, 
                                edge) %>% 
    mutate(agdist = scale(ag_dist),
           agdist_sq = agdist^2,
           agcov = scale(all_ag_pcov_100m),
           agcov_sq = agcov^2,
           aspect = scale(aspect_100m),
           devdist = scale(developed_dist),
           devdist_sq = devdist^2,
           devcov = scale(developed_pcov_1km),
           devcov_sq = devcov^2,
           activity = scale(FPA_TH_activity),
           ripdist = scale(riparian_dist),
           roaddist = scale(road_dist),
           roaddist_sq = roaddist^2,        
           shrubcov = scale(shrub_pcov_500m),
           shrubcov_sq = shrubcov^2,
           slope = scale(slope_1km),
           slope_sq = slope^2,
           RAPtree = scale(RAP_treecover),
           RAPtree_sq = RAPtree^2,
           edge = scale(edge))
  return(dout)
}

# Linear prediction function
modPredict <- function(covs, mod, scale = "exp", keepFirst = FALSE, constrain = NULL){
  mod <- summary(mod)
  coefs <- mod$coefficients[,1]
  se <- mod$coefficients[,2]
  ci <- se*1.96 # ------------------------ adding confidence interval in here in the case we need it later
  if(keepFirst == FALSE){
    coefs <- coefs[-1]
    se <- se[-1]
    ci <- ci[-1]
  }
  
  # Create linear prediction (and transform, as specified in function call)
  pred <- (covs %*% coefs)
  se_up <- (covs %*% (coefs + se))
  se_low <- (covs %*% (coefs - se))
  ci_up <- (covs %*% (coefs + ci))
  ci_low <- (covs %*% (coefs - ci))
  if(!is.null(constrain)) pred <- pConstrain(pred, constrain)
  if(scale == "exp"){
    pred <- exp(pred)
    se_up <- exp(se_up)
    se_low <- exp(se_low)
    ci_up <- exp(ci_up)
    ci_low <- exp(ci_low)
  } 
  
  if(scale == "response"){
    pred <- exp(pred)/(1+exp(pred))
    se_up <- exp(se_up)/(1+exp(se_up))
    se_low <- exp(se_low)/(1+exp(se_low))
    ci_up <- exp(ci_up)/(1+exp(ci_up))
    ci_low <- exp(ci_low)/(1+exp(ci_low))
  } 
  pred_out <- data.frame(pred = pred, se_up = se_up,  se_low = se_low, 
                         ci_up = ci_up, ci_low = ci_low)
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
                       se_col = '#53A583', m_col = 'grey60', scale = "exp", keepFirst = FALSE, constrain = NULL){
  # Prep empty dataframe
  if(keepFirst) d_names <- names(fixef(mod)) else d_names <- names(fixef(mod[-1]))
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


# ------------------------- Run on each age/sex class --------------------------
# Set colors to match SSFs
af_col = '#7f32a8'
am_col = '#e3820b'
df_col = '#d184db'
dm_col = '#e3bf0b'
pooled_col = 'grey50'

# prep datasets
ad_fem <- selected_df %>% filter(age_sex == "Female_adult") %>% datPrep()
ad_male <- selected_df %>% filter(age_sex == "Male_adult") %>% datPrep()
disp_fem <- selected_df %>% filter(age_sex == "Female_disperser") %>% datPrep()
disp_male <- selected_df %>% filter(age_sex == "Male_disperser") %>% datPrep()

# tree
af_tree <- nonLinPlot(mod = ad_fem_models_final_top,
                      cov_name = "RAPtree",
                      orig_name = "RAP_treecover",
                      data = ad_fem,
                      y_name = "Relative probability of selection",
                      keepFirst = TRUE,
                      m_col = af_col,
                      se_col = af_col)

am_tree <- nonLinPlot(mod = ad_male_models_final_top,
                      cov_name = "RAPtree",
                      orig_name = "RAP_treecover",
                      data = ad_male,
                      keepFirst = TRUE,
                      m_col = am_col,
                      se_col = am_col)

df_tree <- nonLinPlot(mod = disp_fem_models_final_top,
                      cov_name = "RAPtree",
                      orig_name = "RAP_treecover",
                      data = disp_fem,
                      x_name = "Annual tree % cover",
                      y_name = "Relative probability of selection",
                      keepFirst = TRUE,
                      m_col = df_col,
                      se_col = df_col)

dm_tree <- nonLinPlot(mod = disp_male_models_final_top,
                      cov_name = "RAPtree",
                      orig_name = "RAP_treecover",
                      data = disp_male,
                      x_name = "Annual tree % cover",
                      keepFirst = TRUE,
                      m_col = dm_col,
                      se_col = dm_col)

pdf(file = 'age-sex-class-RSFs/model-outputs/rsf-tree-nonlin-plot.pdf', height = 6, width = 6)
grid.arrange(af_tree$plot, am_tree$plot, df_tree$plot, dm_tree$plot, nrow = 2)
dev.off()


# shrub
af_shrub <- nonLinPlot(mod = ad_fem_models_final_top,
                       cov_name = "shrubcov",
                       orig_name = "shrub_pcov_500m",
                       data = ad_fem,
                       y_name = "Relative probability of selection",
                       keepFirst = TRUE,
                       x_round = 2,
                       m_col = af_col,
                       se_col = af_col)

df_shrub <- nonLinPlot(mod = disp_fem_models_final_top,
                       cov_name = "shrubcov",
                       orig_name = "shrub_pcov_500m",
                       data = disp_fem,
                       keepFirst = TRUE,
                       x_round = 2,
                       m_col = df_col,
                       se_col = df_col)

dm_shrub <- nonLinPlot(mod = disp_male_models_final_top,
                       cov_name = "shrubcov",
                       orig_name = "shrub_pcov_500m",
                       data = disp_male,
                       x_name = "Shrub % cover",
                       y_name = "Relative probability of selection",
                       keepFirst = TRUE,
                       m_col = dm_col,
                       se_col = dm_col)
              
pdf(file = 'age-sex-class-RSFs/model-outputs/rsf-shrub-nonlin-plot.pdf', height = 6, width = 6)
grid.arrange(af_shrub$plot, df_shrub$plot, dm_shrub$plot, nrow = 2)
dev.off()


# ag dist
af_ag <- nonLinPlot(mod = ad_fem_models_final_top,
                     cov_name = "agdist",
                     orig_name = "ag_dist",
                     data = ad_fem,
                     x_name = "Dist to ag",
                     y_name = "Relative probability of selection",
                     keepFirst = TRUE,
                     x_round = 2,
                     m_col = af_col,
                     se_col = af_col)
pdf(file = 'age-sex-class-RSFs/model-outputs/ag-nonlin-plot.pdf', height = 3, width = 3)
grid.arrange(af_ag$plot, nrow = 1)
dev.off()


# dev dist
af_dev <- nonLinPlot(mod = ad_fem_models_final_top,
                     cov_name = "devdist",
                     orig_name = "developed_dist",
                     data = ad_fem,
                     x_name = "Dist to developed",
                     y_name = "Relative probability of selection",
                     keepFirst = TRUE,
                     x_round = 2,
                     m_col = af_col,
                     se_col = af_col)
pdf(file = 'age-sex-class-RSFs/model-outputs/rsf-dev-nonlin-plot.pdf', height = 3, width = 3)
grid.arrange(af_dev$plot, nrow = 1)
dev.off()


# slope
af_slope <- nonLinPlot(mod = ad_fem_models_final_top,
                       cov_name = "slope",
                       orig_name = "slope_1km",
                       data = ad_fem,
                       y_name = "Relative probability of selection",
                       keepFirst = TRUE,
                       x_round = 0,
                       m_col = af_col,
                       se_col = af_col)

am_slope <- nonLinPlot(mod = ad_male_models_final_top,
                       cov_name = "slope",
                       orig_name = "slope_1km",
                       x_name = "Slope",
                       data = ad_male,
                       keepFirst = TRUE,
                       x_round = 0,
                       m_col = am_col,
                       se_col = am_col)

dm_slope <- nonLinPlot(mod = disp_male_models_final_top,
                       cov_name = "slope",
                       orig_name = "slope_1km",
                       data = disp_male,
                       keepFirst = TRUE,
                       x_round = 0,
                       m_col = dm_col,
                       se_col = dm_col)
pdf(file = 'age-sex-class-RSFs/model-outputs/rsf-slope-nonlin-plot.pdf', height = 3, width = 9)
grid.arrange(af_slope$plot, am_slope$plot, dm_slope$plot, nrow = 1)
dev.off()