## ---------------------------
##
## Script name: age-sex-pooled-05-RSF-vis.R
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

# Load final pooled model
load(file = "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/top-final-age-sex-pooled-RSF.rda", verbose = TRUE)

# Prep model coefficients
pooledcoef <- broom.mixed::tidy(pooled_age_sex_models_final_top) %>% filter(effect == 'fixed')
pooledcoef <- transform(pooledcoef,
                        term=sprintf("%s.%s", effect, term))[-1,]

# save CSVs for later
write_csv(pooledcoef, "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/age-sex-pooled-RSF-coefficients.csv")

# ------------------------------------ Plot ------------------------------------
pdf(file = "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/age-sex-pooled-rsf-coefs.pdf", width = 5, height = 4)
pooled_coef_plot <- dwplot(pooledcoef,
                           ci = 0.95,
                           dot_args = list(size = 2.5, shape = 22, fill = 'grey50', color = 'black'),
                           whisker_args = list(size = 1, color = "black"),
                           vline = geom_vline(xintercept = 0,
                                              colour = "grey60",
                                              linetype = 2)) %>% 
  relabel_predictors(fixed.edge = "Annual edge % cover",
                     fixed.RAPtree = "Annual tree % cover",
                     fixed.RAPtree_sq = "Annual tree^2",
                     fixed.shrubcov = "Shrub % cover",
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
pooled_coef_plot
dev.off()


# ------------------------------------------------------------------------------ 
# ---------------------- Plot non-linear covariate effects ---------------------
# ------------------------------------------------------------------------------ 

# Load raw data table
load("age-sex-class-RSFs/age-sex-pooled-RSF/age_sex_pooled_RSF_covariate_extract_all.rda")

# --------------------------------- Functions ---------------------------------

# Prep pooled dataset
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
  ci <- se*1.96 # ------------------------ adding confience interval in here in the case we need it later
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
                       se_col = '#FF7373', m_col = 'grey60', scale = "exp", keepFirst = FALSE, constrain = NULL){
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


# ---------------------------- Run on pooled model -----------------------------

# prep dataset
pooled <- selected_df %>% datPrep()

# tree
p_tree <- nonLinPlot(mod = pooled_age_sex_models_final_top,
                     cov_name = "RAPtree",
                     orig_name = "RAP_treecover",
                     data = pooled,
                     y_name = "Relative probability of selection",
                     x_name = "Annual tree % cover",
                     keepFirst = TRUE) 
pdf(file = "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/pooled-rsf-tree-nonlin.pdf", width = 5, height = 5)
p_tree$plot
dev.off()

# ag dist
p_agdist <- nonLinPlot(mod = pooled_age_sex_models_final_top,
                   cov_name = "agdist",
                   orig_name = "ag_dist",
                   data = pooled,
                   x_div = 1000,
                   # x_round = 2,
                   x_name = "Distance to agriculture (km)",
                   y_name = "Relative probability of selection",
                   keepFirst = TRUE)
pdf(file = "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/pooled-rsf-agdist-nonlin.pdf", width = 5, height = 5)
p_agdist$plot
dev.off()

# dev dist
p_devdist <- nonLinPlot(mod = pooled_age_sex_models_final_top,
                    cov_name = "devdist",
                    orig_name = "developed_dist",
                    data = pooled,
                    x_div = 1000,
                    # x_round = 2,
                    x_name = "Distance to development (km)",
                    y_name = "Relative probability of selection",
                    keepFirst = TRUE)
pdf(file = "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/pooled-rsf-devdist-nonlin.pdf", width = 5, height = 5)
p_devdist$plot
dev.off()

# road dist
p_roaddist <- nonLinPlot(mod = pooled_age_sex_models_final_top,
                     cov_name = "roaddist",
                     orig_name = "road_dist",
                     data = pooled,
                     x_div = 1000,
                     # x_round = 2,
                     x_name = "Distance to road (km)",
                     y_name = "Relative probability of selection",
                     keepFirst = TRUE)
pdf(file = "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/pooled-rsf-roaddist-nonlin.pdf", width = 5, height = 5)
p_roaddist$plot
dev.off()

# slope dist
p_slope <- nonLinPlot(mod = pooled_age_sex_models_final_top,
                      cov_name = "slope",
                      orig_name = "slope_1km",
                      data = pooled,
                      # x_round = 2,
                      x_name = "Slope (degrees)",
                      y_name = "Relative probability of selection",
                      keepFirst = TRUE)
pdf(file = "age-sex-class-RSFs/age-sex-pooled-RSF/model-outputs/pooled-rsf-slope-nonlin.pdf", width = 5, height = 5)
p_slope$plot
dev.off()

