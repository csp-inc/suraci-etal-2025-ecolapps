## ---------------------------
##
## Script name: sighting-conflict-analysis-glm.R
##
## Purpose of script: A simple modeling exercise to evaluate the explanatory effect of age/sex-class specific movement SSFs and kill-site RSFs on observed instances of reported puma-human conflict
##
## Author: Patrick Freeman
##
## Date Created: 10/23/2023
## Date last updated: 03/06/2024
##
## Email contact: patrick[at]csp-inc.org
##
## ---------------------------
##
## Notes: 
##  

library(sf)
library(terra)
library(exactextractr)
library(predicts)
library(tidyterra)
library(tidyverse)
library(ggthemes)
library(tidymodels)
library(broom)
library(purrr)
library(dotwhisker)
library(lubridate)
library(corrr)
library(ggeffects)
library(ggspatial)
library(ResourceSelection)
library(maptiles)


# Utility Functions -------------------------------------------------------
### Rescale rasters to 0-1
rescale01 <- function(r){
  mm <- minmax(r)
  rout <- (r - mm[1])/(mm[2] - mm[1])
  return(rout)
}

### Prepare model data
### Supply with underlying raster surface, the sampling domain polygon, the conflict data set, and a model name 
### Prepares a ready-to-go dataframe for a binomial GLM 
prepare_model_data <- function(raster=raster,
                               pt_data = conflict_background,
                               model_name = "model_name"){
  
  data_prep_sample <- terra::extract(raster, pt_data, bind=T, fun="mean", na.rm=T) %>%
    as_tibble() %>%
    dplyr::mutate(point_type = case_when(point_type == "conflict" ~ 1,
                                         point_type == "background" ~ 0)) %>%
    rename(conflict = point_type) %>%
    #dplyr::rename(predictor = model_name) %>%
    dplyr::mutate(id = as.character(row_number(.))) %>%
    dplyr::mutate(join_col = paste(id, conflict, ReportDate, month, sep="_")) %>%
    dplyr::select(-c(conflict,ReportDate,month, id))
  
  return(data_prep_sample)
  
}

# Load conflict data ------------------------------------------------------

conflict_raw <- st_read("data/cougage-ALL-conflict-2016to2022/cougar_ALL_conflict_2016to2022.shp") %>%
  st_transform(., crs=26910)

conflict_sight_confront <- conflict_raw %>%
  dplyr::filter(IncidentTy %in% c("Sighting", "Confrontation"))

# Load raster data ------------------------------------------------------

all_cats_ssf <- rast("output/all-cats-ssf-prediction-30m.tif")
names(all_cats_ssf) <- "all_cats_SSF"

ad_m_ssf <- rast("output/ad-male-ssf-prediction-30m.tif")
names(ad_m_ssf) <- "ad_m_SSF"

disp_m_ssf <- rast("output/disp-male-ssf-prediction-30m.tif")
names(disp_m_ssf) <- "disp_m_SSF"

disp_f_ssf <- rast("output/disp-fem-ssf-prediction-30m.tif")
names(disp_f_ssf) <- "disp_f_SSF"

ad_f_ssf <- rast("output/ad-fem-ssf-prediction-30m.tif")
names(ad_f_ssf) <- "ad_f_SSF"


all_cats_rsf <- rast("output/age-sex-pooled-exp-aoi-30.tif")
names(all_cats_rsf) <- "all_cats_RSF"

ad_m_rsf <- rast("output/ad-male-exp-aoi-30.tif")
names(ad_m_rsf) <- "ad_m_RSF"

disp_m_rsf <- rast("output/disp-male-exp-aoi-30.tif")
names(disp_m_rsf) <- "disp_m_RSF"

disp_f_rsf <- rast("output/disp-fem-exp-aoi-30.tif")
names(disp_f_rsf) <- "disp_f_RSF"

ad_f_rsf <- rast("output/ad-fem-exp-aoi-30.tif")
names(ad_f_rsf) <- "ad_f_RSF"


# Load sampling domain ----------------------------------------------------

sampling_domain <- st_read("data/conflict-sampling-domain/all_ag_VIIRS_dev_full_yr_comb_buffer3km_clipland.shp")

# Tidy livestock conflict database ------------------------------------------------------
sighting_conflict_tidy <- conflict_sight_confront %>%
  dplyr::select(geometry, IncidentTy, ReportDate) %>%
  dplyr::rename(point_type = IncidentTy) %>% 
  dplyr::relocate(point_type, .before=everything()) %>%
  dplyr::mutate(month = month(ReportDate)) %>%
  dplyr::mutate(point_type = "conflict")


# Create 'background' sample of points using ------------------------------------------------------
#### Using the backgroundSample function from the 'predicts' package

conflict_n <- nrow(sighting_conflict_tidy )

set.seed(1234)
### Create mask to remove possibiilty of sampling in the water
mask <- terra::mask(ad_m_ssf, sampling_domain)
### Create background samples 
background_pts <- backgroundSample(mask=mask, n=10*conflict_n, p=vect(sighting_conflict_tidy), excludep=T, tryf=9) %>%
  as_tibble() %>%
  st_as_sf(coords=c(c("x", "y")), crs=crs(sighting_conflict_tidy)) %>%
  dplyr::mutate(point_type="background") %>%
  dplyr::relocate(point_type, .before=everything())

# Create final conflict-background dataset ------------------------------------------------------
all_conflict_bg <- bind_rows(sighting_conflict_tidy, background_pts) %>%
  st_buffer(., 100)


# Plot conflict and background points for visual check  ------------------------------------------------------
op <- get_tiles(sampling_domain, crop=T, zoom=9, provider="Esri.WorldImagery")

ggplot() + 
  geom_spatraster_rgb(data = op, maxcell = Inf) + 
  scale_fill_viridis_c() + 
  geom_spatvector(data=vect(sampling_domain), fill=NA, color="gold", lwd=0.5) + 
  geom_sf(data=st_as_sf(all_conflict_bg), aes(color=point_type), size=0.75, show.legend="point") + 
  scale_color_manual(values=c("darkgrey", "red")) + 
  scale_shape_manual(values = c(16, 16)) + 
  guides(color = guide_legend(override.aes = list(size = 10))) + 
  annotation_scale() + 
  labs(color="Point type",
       title="Livestock conflict and background points") + 
  theme(legend.position="right")

# Prepare Model Data  ------------------------------------------------------
### Create list of SSF rasters
puma_ssf_list <- list(all_cats_ssf, ad_m_ssf, disp_m_ssf, ad_f_ssf, disp_f_ssf)
names(puma_ssf_list) <- c("all_cats_ssf", "ad_m_SSF", "disp_m_SSF", "ad_f_SSF", "disp_f_SSF")
### Create a list of model names from the named list
ssf_model_names <- names(puma_ssf_list)

### Clip to sampling domain
samp_domain_list <- list(sampling_domain, sampling_domain, sampling_domain, sampling_domain, sampling_domain)


puma_ssf_mask_list <- map2(puma_ssf_list, samp_domain_list, ~terra::mask(.x, .y))


### Create list of kill site RSF rasters
puma_kill_rsf_list <- list(all_cats_rsf, ad_m_rsf, disp_m_rsf, ad_f_rsf, disp_f_rsf)
names(puma_kill_rsf_list) <- c("all_cats_rsf", "ad_m_rsf", "disp_m_rsf", "ad_f_rsf", "disp_f_rsf")
rsf_model_names <- names(puma_kill_rsf_list)

### Mask kill site RSF rasters 
puma_rsf_mask_list <- map2(puma_kill_rsf_list, samp_domain_list, ~terra::mask(.x, .y))

### Prepare the model data for all surfaces
puma_ssf_stack <- rast(puma_ssf_mask_list)
ssf_data_extract <- terra::extract(puma_ssf_stack, all_conflict_bg, bind=T, fun="mean", na.rm=T) %>%
  as_tibble() %>%
  dplyr::mutate(point_type = case_when(point_type == "conflict" ~ 1,
                                       point_type == "background" ~ 0)) %>%
  rename(conflict = point_type) %>%
  dplyr::mutate(id = as.character(rownames(.)))


puma_rsf_stack <- rast(puma_rsf_mask_list)
rsf_data_extract <- terra::extract(puma_rsf_stack, all_conflict_bg, bind=T, fun="mean", na.rm=T) %>%
  as_tibble() %>%
  dplyr::mutate(point_type = case_when(point_type == "conflict" ~ 1,
                                       point_type == "background" ~ 0)) %>%
  rename(conflict = point_type) %>%
  dplyr::mutate(id = as.character(rownames(.))) 


# Join and Z-scale Prepared Model Data  -----------------------------------------------

sighting_conflict_prepared_data <- left_join(ssf_data_extract, rsf_data_extract, by=c("id", "conflict", "ReportDate","month")) %>%
  dplyr::relocate(id, .before=everything()) %>%
  mutate_at(c(5:14), funs(c(scale(.))))


# Create model datasets  --------------------------------------------------

ad_m_sighting_conflict <- sighting_conflict_prepared_data  %>% 
  dplyr::select(id, conflict, ReportDate, month, ad_m_SSF, ad_m_rsf) 

disp_m_sighting_conflict <- sighting_conflict_prepared_data  %>% 
  dplyr::select(id, conflict, ReportDate, month, disp_m_SSF, disp_m_rsf)

ad_f_sighting_conflict <- sighting_conflict_prepared_data  %>% 
  dplyr::select(id, conflict, ReportDate, month, ad_f_SSF, ad_f_rsf)

disp_f_sighting_conflict <- sighting_conflict_prepared_data  %>% 
  dplyr::select(id, conflict, ReportDate, month, disp_f_SSF, disp_f_rsf)

allcats_sighting_conflict <- sighting_conflict_prepared_data  %>%
  dplyr::select(id, conflict, ReportDate, month, all_cats_ssf, all_cats_rsf)

# Build GLMs for each age-sex class - combining the RSF and SSF for each combination into a single model with two predictors ------------------------------------------------------

ad_m_sighting_model <- glm(conflict ~ ad_m_SSF + ad_m_rsf, data=ad_m_sighting_conflict, family="binomial")
disp_m_sighting_model <- glm(conflict ~ disp_m_SSF + disp_m_rsf, data=disp_m_sighting_conflict, family="binomial")
ad_f_sighting_model <- glm(conflict ~ ad_f_SSF + ad_f_rsf, data=ad_f_sighting_conflict, family="binomial")
disp_f_sighting_model <- glm(conflict ~ disp_f_SSF + disp_f_rsf, data=disp_f_sighting_conflict, family="binomial")
all_cats_sighting_model <- glm(conflict ~ all_cats_ssf + all_cats_rsf, data= allcats_sighting_conflict, family="binomial")



# Tidy up GLM outputs ------------------------------------------------------

### 'Tidy output with coefficient estimates'
### Age/sex class Models 
ad_m_sighting_tidy_out <- tidy(ad_m_sighting_model, conf.int = TRUE, conf.level = 0.95, exponentiate=F) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(varType = case_when(str_detect(term, "SSF") == T ~ "ssf",
                                    str_detect(term, "rsf") == T ~ "rsf")) %>% 
  dplyr::mutate(ASClass = "adult_male")

disp_m_sighting_tidy_out <- tidy(disp_m_sighting_model, conf.int = TRUE, conf.level = 0.95, exponentiate=F) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>% 
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(varType = case_when(str_detect(term, "SSF") == T ~ "ssf",
                                    str_detect(term, "rsf") == T ~ "rsf")) %>% 
  dplyr::mutate(ASClass = "disp_male")

ad_f_sighting_tidy_out <- tidy(ad_f_sighting_model, conf.int = TRUE, conf.level = 0.95, exponentiate=F) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>% 
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(varType = case_when(str_detect(term, "SSF") == T ~ "ssf",
                                    str_detect(term, "rsf") == T ~ "rsf")) %>%
  dplyr::mutate(ASClass = "adult_female")

disp_f_sighting_tidy_out <- tidy(disp_f_sighting_model, conf.int = TRUE, conf.level = 0.95, exponentiate=F) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>% 
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(varType = case_when(str_detect(term, "SSF") == T ~ "ssf",
                                    str_detect(term, "rsf") == T ~ "rsf")) %>% 
  dplyr::mutate(ASClass = "disp_female")

allcats_sighting_tidy_out <- tidy(all_cats_sighting_model, conf.int = TRUE, conf.level = 0.95, exponentiate=F) %>%
  dplyr::select(term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>% 
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(varType = case_when(str_detect(term, "ssf") == T ~ "ssf",
                                    str_detect(term, "rsf") == T ~ "rsf")) %>% 
  dplyr::mutate(ASClass = "all_cats")

# Bind GLM Coefficient Estimates into single table ------------------------------------------------------

all_sighting_model_coefs <- bind_rows(list(ad_m_sighting_tidy_out, disp_m_sighting_tidy_out, ad_f_sighting_tidy_out, disp_f_sighting_tidy_out, allcats_sighting_tidy_out))

all_sighting_model_coefs $ASClass <- factor(all_sighting_model_coefs$ASClass, levels=rev(c("all_cats", "adult_female", "adult_male", "disp_female", "disp_male")), ordered=T)


all_sighting_model_coefs <- all_sighting_model_coefs %>%
  dplyr::mutate(
    size_multiplier = case_when(ASClass == "all_cats" ~ 1.05,
                                ASClass != "all_cats" ~ 1))

# Plot model coefficients ------------------------------------------------------
# Set colors
af_col = '#7f32a8'
am_col = '#e3820b'
df_col = '#d184db'
dm_col = '#e3bf0b'
all_col = 'grey'

### Make forest plot 
sightingdw_plot <- dwplot(arrange(all_sighting_model_coefs, desc(ASClass)),
       ci = 0.95,
       dot_args = list(size = 6, 
                       aes(color = ASClass,
                           shape=varType)),
       whisker_args = list(size = 1, 
                           color = "black"),
       vline = geom_vline(xintercept = 0,
                          colour = "grey60",
                          linetype = 2)
) +
  xlim(-2,2) + 
  scale_color_manual(values=c(all_col, am_col, dm_col, af_col, df_col), breaks=c("all_cats", "adult_male", "disp_male", "adult_female", "disp_female"), labels=c("all cats", "adult male", "dispersing male", "adult female", "dispersing female"), name="Age/Sex Class") + 
  scale_shape_discrete(name  ="Predictor Type", breaks = c("rsf", "ssf"), labels=c("Kill site RSF", "Movement SSF")) + # breaks assign shapes
  theme_clean() + 
  labs(x="Coefficient Estimate", y="Puma Model", title="Sighting and Confrontation Conflicts Only") + 
  theme(text=element_text(size=24, face="bold"),
        axis.text.y = element_blank())


(sighting_conflict_coef_plot <- ggplot(data = all_sighting_model_coefs, 
                                        aes(x = ASClass, 
                                            y = estimate, 
                                            fill = ASClass, 
                                            shape=varType),
                                        color="black") + 
    geom_pointrange(aes(ymax = conf.high, 
                        ymin = conf.low,
                        size= size_multiplier)) +
    scale_size(range = c(1, 1.25)) + 
    scale_fill_manual(values=c(all_col, af_col, am_col, df_col, dm_col), 
                      breaks=c("all_cats", "adult_female", "adult_male", "disp_female", "disp_male"), 
                      labels=c("pooled", "adult female", "adult male", "dispersing female", "dispersing male"),
                      name="Age/Sex Class") + 
    scale_shape_manual(name  ="Model Type", 
                       breaks = c("rsf", "ssf"), # breaks assign shapes
                       values = c(22, 21),
                       labels=c("Kill site RSF", "Movement SSF")) + 
    coord_flip() + 
    geom_hline(yintercept = 0,
               colour = "grey60",
               linetype = 2) + 
    theme_bw() + 
    guides(fill = guide_legend(override.aes = list(shape = 21, size=1.5) ),
           shape = guide_legend(override.aes = list(fill = "black", size=1.5)) ,
           size="none"
    ) + 
    theme(
      text=element_text(size=14, face="bold"),
      axis.text.y = element_blank(),
      axis.title.y = element_blank()))

pdf(file = '00-out/sighting-only-conflict-coefs-20240306.pdf', height = 6, width = 6)
grid.arrange(sighting_conflict_coef_plot, nrow = 1)
dev.off()
