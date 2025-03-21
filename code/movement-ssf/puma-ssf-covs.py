# -*- coding: utf-8 -*-
import ee
import math

# Initialize Earth Engine
ee.Initialize()

# Define AOI
# AOI = ee.FeatureCollection('projects/GEE_CSP/Panthera/puma_AOI')
AOI = ee.FeatureCollection('projects/GEE_CSP/Panthera/refined_AOI_03Feb2023')

# Export scale and production
scale = 30
projection = ee.Projection('ESPG:26910') # NAD83 UTM Zone 10N

# Choose radii for summarizing covariates (three alternatives to test in model)
rad_large = 1000 # need to confirm these radii
name_large = "_1km"
# rad_med = 500 # need to confirm these radii 
# name_med = "_500m"
rad_small = 150 # need to confirm these radii
name_small = "_150m"

#---------------------------------------------------------------
# FUNCTIONS

# Focal mean
def focal_mean(image, radius, unit, name):
    names = image.bandNames().getInfo()
    new_names = [s + name for s in names]
    return image.reduceNeighborhood(kernel = ee.Kernel.circle(radius, unit),
                                    reducer = ee.Reducer.mean()).rename(new_names)

# Focal sum
def focal_sum(image, radius, unit):
    kernel = ee.Kernel.circle(radius, unit, False)
    reducer = ee.Reducer.sum()
    result = image.reduceNeighborhood(reducer,kernel)
    return result

# Focal count
def focal_count(image, radius, unit):
    return image.reduceNeighborhood(kernel = ee.Kernel.circle(radius, unit, False),
                                    reducer = ee.Reducer.count())

# Percent cover
def percent_cov(image, radius, unit, name):
    names = image.bandNames().getInfo()
    new_names = [s + name for s in names]
    isum = focal_sum(image, radius, unit)
    icount = focal_count(image, radius, unit)
    return isum.divide(icount).rename(new_names)

# Vector ruggedness measure
def compute_vrm(slope_img, aspect_img, radius, units):
    slope_sine = slope_img.sin()
    x_sum_sq = focal_sum(slope_sine.multiply(aspect_img.sin()), radius, units).pow(2)
    y_sum_sq = focal_sum(slope_sine.multiply(aspect_img.cos()), radius, units).pow(2)
    z_sum_sq = focal_sum(slope_img.cos(), radius, units).pow(2)
    n = focal_sum(ee.Image(1), radius, units)
    r = x_sum_sq.add(y_sum_sq).add(z_sum_sq).sqrt()
    vrm_img = ee.Image(1).subtract(r.divide(n))
    return vrm_img

# Topographic position index
def topo_position(elev_img, radius, unit):
    En = elev_img.focalMean(radius, 'circle', unit)
    tpi_out = elev_img.subtract(En).abs()
    return tpi_out

# Iterate over all study years for seasonal EVI
def springEVI(year):
    startDate = ee.Date.fromYMD(year,2,15)
    endDate = ee.Date.fromYMD(year,5,15)
    composite_i = EVIimageCollection.filterDate(startDate, endDate).median()
    return composite_i

def summerEVI(year):
    startDate = ee.Date.fromYMD(year,5,15)
    endDate = ee.Date.fromYMD(year,8,15)
    composite_i = EVIimageCollection.filterDate(startDate, endDate).median()
    return composite_i

def fallEVI(year):
    startDate = ee.Date.fromYMD(year,8,15)
    endDate = ee.Date.fromYMD(year,11,15)
    composite_i = EVIimageCollection.filterDate(startDate, endDate).median()
    return composite_i

def winterEVI(year):
    startDate = ee.Date.fromYMD(year,11,15)
    endDate = ee.Date.fromYMD(ee.Number(year).add(1),2,15)
    composite_i = EVIimageCollection.filterDate(startDate, endDate).median()
    return composite_i

#---------------------------------------------------------------
# TERRAIN VARIABLES
#---------------------------------------------------------------

# Create some topographic viables
dsm = ee.ImageCollection("JAXA/ALOS/AW3D30/V3_2") # Digital surface model. Native res 30m
projElev = dsm.first().select(0).projection()
elev = dsm.select("DSM").mosaic().setDefaultProjection(projElev).rename('elevation').clip(AOI) # Elevation
slope = ee.Terrain.slope(elev).rename('slope').clip(AOI) # Slope
aspect = ee.Terrain.aspect(elev).rename('aspect').clip(AOI) # Aspect

# Create vector ruggedness metric
window_radius = 1000
slopeRad = slope.multiply(ee.Number(math.pi).divide(180))
aspectRad = aspect.multiply(ee.Number(math.pi).divide(180))
vrm = compute_vrm(slopeRad, aspectRad, window_radius, "meters").rename('vrm').clip(AOI)
# Create TPI
tpi = topo_position(elev, window_radius, "meters").rename("tpi").clip(AOI)

# Combine and smooth
terrain_all = ee.Image([elev, slope, aspect, vrm, tpi])
terrain_large = focal_mean(terrain_all, rad_large, "meters", name_large)
terrain_small = focal_mean(terrain_all, rad_small, "meters", name_small)

#---------------------------------------------------------------
# ENVIRO COVARIATES
#---------------------------------------------------------------

# setting year lists for seasonal enviro covariates
stepList = ee.List.sequence(2018,2022) # Make a list of years (based on when kill site data available), then for each year filter the collection for seasonal time frames
stepListWinter = ee.List.sequence(2018,2021)

# EVI
EVIimageCollection = ee.ImageCollection('MODIS/MCD43A4_006_EVI').filterBounds(AOI) # Import MODIS EVI
springFilterCollection = stepList.map(springEVI)
summerFilterCollection = stepList.map(summerEVI)
fallFilterCollection = stepList.map(fallEVI)
winterFilterCollection = stepListWinter.map(winterEVI)
springEVIComposites = ee.ImageCollection(springFilterCollection).select('EVI').mosaic().setDefaultProjection(projElev).clip(AOI).rename('spring_EVI')
summerEVIComposites = ee.ImageCollection(summerFilterCollection).select('EVI').mosaic().setDefaultProjection(projElev).clip(AOI).rename('summer_EVI')
fallEVIComposites = ee.ImageCollection(fallFilterCollection).select('EVI').mosaic().setDefaultProjection(projElev).clip(AOI).rename('fall_EVI')
winterEVIComposites = ee.ImageCollection(winterFilterCollection).select('EVI').mosaic().setDefaultProjection(projElev).clip(AOI).rename('winter_EVI')

## smoothing EVI composites
EVI_all = ee.Image([springEVIComposites, summerEVIComposites, fallEVIComposites, winterEVIComposites])
EVI_large = focal_mean(EVI_all, rad_large, "meters", name_large)
EVI_small = focal_mean(EVI_all, rad_small, "meters", name_small)

#---------------------------------------------------------------
# LAND COVER
#---------------------------------------------------------------
# PREP DATA
lc = ee.ImageCollection('USGS/NLCD_RELEASES/2019_REL/NLCD') # Import NLCD
nlcd2019 = lc.filter(ee.Filter.eq('system:index', '2019')).first() # Filter the collection to the 2019 product
landcover = nlcd2019.select('landcover').setDefaultProjection(projElev).clip(AOI) # Select land cover band
forest = ee.Image(0).where(landcover.gte(41).And(landcover.lte(43)), 1).rename('forest').setDefaultProjection(projElev).clip(AOI)
shrub = landcover.eq(52).rename('shrub').setDefaultProjection(projElev).clip(AOI)
grass = landcover.eq(71).rename('grass').setDefaultProjection(projElev).clip(AOI)
impervious = nlcd2019.select('impervious').setDefaultProjection(projElev).clip(AOI).rename('impervious') 

# Distance to NHD ripariran zones
riparian = ee.Image('projects/GEE_CSP/Panthera/EPA-riparian-zones-2019') # Import riparian zone image

# Agriculture - adding NLCD crop and pasture layers with AFT FUT crop layer
## NLCD crop + pasture
nlcd_ag = ee.Image(0).where(landcover.gte(81).And(landcover.lte(82)), 1).rename('nlcd_ag').setDefaultProjection(projElev).clip(AOI)
## AFT FUT crop
aft_lcov = ee.Image('projects/GEE_CSP/AFT_FUT/LCU_2015_v6')
aft_ag = aft_lcov.expression(
    "(b(0) == 1101 | b(0) == 1) ? 1" +
      ": (b(0) == 1102 | b(0) == 2) ? 2" +
        ": (b(0) == 1103 | b(0) == 3) ? 3" +
          ": (b(0) == 1105 | b(0) == 5) ? 5" +
            ": 0")
# Binary layers for ag/not ag
isCrop = aft_ag.eq(1).rename('aft_crop')
# Combine binary layers
all_ag = nlcd_ag.add(isCrop).rename('ag').gte(1)

# SUMMARIZE
# Focal
lc_all = ee.Image([forest, shrub, grass, impervious, all_ag])
lc_large = percent_cov(lc_all, rad_large, 'meters', '_pcov' + name_large)
lc_small = percent_cov(lc_all, rad_small, 'meters', '_pcov' + name_small)

# Distance
rip_dist = riparian.distance(ee.Kernel.euclidean(50000, 'meters'), skipMasked = False).clip(AOI).rename('riparian_dist')
ag_dist = all_ag.distance(ee.Kernel.euclidean(50000, 'meters'), skipMasked = False).clip(AOI).rename('ag_dist')

#---------------------------------------------------------------
# FORESTS/FORESTRY
#---------------------------------------------------------------
# COLLECT DATA

# RAP % tree cover
RAP = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3').toBands()
RAP2017 = RAP.select('2017_TRE').setDefaultProjection(projElev).clip(AOI).rename('2017_trCover')
RAP2018 = RAP.select('2018_TRE').setDefaultProjection(projElev).clip(AOI).rename('2018_trCover')
RAP2019 = RAP.select('2019_TRE').setDefaultProjection(projElev).clip(AOI).rename('2019_trCover')
RAP2020 = RAP.select('2020_TRE').setDefaultProjection(projElev).clip(AOI).rename('2020_trCover')
RAP2021 = RAP.select('2021_TRE').setDefaultProjection(projElev).clip(AOI).rename('2021_trCover')
RAP2022 = RAP.select('2022_TRE').setDefaultProjection(projElev).clip(AOI).rename('2022_trCover')

RAP_forest2017 = RAP2017.updateMask(RAP2017.gte(ee.Number(45))).unmask().gt(0).rename('RAP2017_forest45pcnt')
RAP_forest2018 = RAP2018.updateMask(RAP2018.gte(ee.Number(45))).unmask().gt(0).rename('RAP2018_forest45pcnt')
RAP_forest2019 = RAP2019.updateMask(RAP2019.gte(ee.Number(45))).unmask().gt(0).rename('RAP2019_forest45pcnt')
RAP_forest2020 = RAP2020.updateMask(RAP2020.gte(ee.Number(45))).unmask().gt(0).rename('RAP2020_forest45pcnt')
RAP_forest2021 = RAP2021.updateMask(RAP2021.gte(ee.Number(45))).unmask().gt(0).rename('RAP2021_forest45pcnt')
RAP_forest2022 = RAP2022.updateMask(RAP2022.gte(ee.Number(45))).unmask().gt(0).rename('RAP2022_forest45pcnt')

## Forest edge (created externally???)
edge_2017 = ee.Image('projects/GEE_CSP/Panthera/2017_TRE_edge').setDefaultProjection(projElev).clip(AOI).rename('2017edge')
edge_2018 = ee.Image('projects/GEE_CSP/Panthera/2018_TRE_edge').setDefaultProjection(projElev).clip(AOI).rename('2018edge')
edge_2019 = ee.Image('projects/GEE_CSP/Panthera/2019_TRE_edge').setDefaultProjection(projElev).clip(AOI).rename('2019edge')
edge_2020 = ee.Image('projects/GEE_CSP/Panthera/2020_TRE_edge').setDefaultProjection(projElev).clip(AOI).rename('2020edge')
edge_2021 = ee.Image('projects/GEE_CSP/Panthera/2021_TRE_edge').setDefaultProjection(projElev).clip(AOI).rename('2021edge')
edge_2022 = ee.Image('projects/GEE_CSP/Panthera/2022_TRE_edge').setDefaultProjection(projElev).clip(AOI).rename('2022edge')

## Disturbed edge dataset (created externally??)
disturbed_edge = ee.Image('projects/GEE_CSP/Panthera/FPA_TH_disturbed_invert_edge_forestmask_clean').setDefaultProjection(projElev).clip(AOI).rename('disturbed_edge')

# Combine and summarize
forest_all = ee.Image([disturbed_edge,RAP2017,RAP2018,RAP2019,RAP2020,RAP2021,RAP2022,
                       edge_2017,edge_2018,edge_2019,edge_2020,edge_2021,edge_2022])
forest_large = percent_cov(forest_all, rad_large, 'meters', '_pcov' + name_large)
forest_small = percent_cov(forest_all, rad_small, 'meters', '_pcov' + name_small)


# Forestry: Forest Practice Applications and Timber Harvest
FPA_TH_activity = ee.Image('projects/GEE_CSP/Panthera/FPA_TH_polys_active_forestry_final').setDefaultProjection(projElev).clip(AOI).rename('forestry')
forestry_large = percent_cov(FPA_TH_activity, rad_large, 'meters', '_pcov' + name_large)
forestry_small = percent_cov(FPA_TH_activity, rad_small, 'meters', '_pcov' + name_small)

#---------------------------------------------------------------
# HUMAN INFLUENCE
#---------------------------------------------------------------
# COLLECT DATASETS
# Human population density
hp = ee.ImageCollection("CIESIN/GPWv411/GPW_Population_Density").first()
human_pop = hp.select('population_density').setDefaultProjection(projElev).clip(AOI).rename('human_pop_den')
# Artificial night light
viirs = ee.ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG").filter(ee.Filter.date('2021-01-01', '2021-05-31'))
nighttime = viirs.select('avg_rad').mosaic().setDefaultProjection(projElev).clip(AOI).rename('nightlight')
# NLCD Developed
low_med_high_dev = ee.Image(0).where(landcover.gte(22).And(landcover.lte(24)), 1).rename('developed').setDefaultProjection(projElev).clip(AOI)
# Roads
roads = ee.FeatureCollection('projects/GEE_CSP/Panthera/WSDOT-Functional_Class_Data_for_State_Routes')
rd_img = roads.reduceToImage(properties = ['FederalFun'], reducer = ee.Reducer.max()).clip(AOI)

#SUMMARIZE
# Focal
human_all = ee.Image([human_pop, nighttime, low_med_high_dev])
human_large = focal_mean(human_all, rad_large, "meters", name_large)
human_small = focal_mean(human_all, rad_small, "meters", name_small)

# Distance
rd_dist = rd_img.distance(ee.Kernel.euclidean(50000, 'meters'), skipMasked = False).clip(AOI).rename('road_dist')
dev_dist = low_med_high_dev.distance(ee.Kernel.euclidean(50000, 'meters'), skipMasked = False).clip(AOI).rename('developed_dist')

#---------------------------------------------------------------
# COMBINE AND EXPORT
#---------------------------------------------------------------
focal_covs = ee.Image([terrain_all, terrain_small, terrain_large,
                       EVI_all, EVI_small, EVI_large,
                       lc_all, lc_small, lc_large,
                       forest_all, forest_small, forest_large, FPA_TH_activity,
                       human_all, human_small, human_large]).float().clip(AOI)

dist_covs = ee.Image([rip_dist, ag_dist, rd_dist, dev_dist]).float().clip(AOI)

forestry_covs = ee.Image([forestry_large, forestry_small]).float().clip(AOI)

# Export focal covs
task1 = ee.batch.Export.image.toCloudStorage(image = focal_covs,
                                             bucket = 'puma-ssf-analysis',
                                             description = 'puma-ssf-focal-covs-' + str(scale) + "m",
                                             scale = scale,
                                             region = AOI.geometry(),
                                             maxPixels = 1e13,
                                             crs = "EPSG:26910")
task1.start()


## Export dist covs
task2 = ee.batch.Export.image.toCloudStorage(image = dist_covs,
                                             bucket = 'puma-ssf-analysis',
                                             description = 'puma-ssf-dist-covs-270m',
                                             scale = 270,
                                             region = AOI.geometry(),
                                             maxPixels = 1e13,
                                             crs = "EPSG:26910")
task2.start()

# Export focal covs
task3 = ee.batch.Export.image.toCloudStorage(image = forestry_covs,
                                             bucket = 'puma-ssf-analysis',
                                             description = 'puma-ssf-forestry-covs-' + str(scale) + "m",
                                             scale = scale,
                                             region = AOI.geometry(),
                                             maxPixels = 1e13,
                                             crs = "EPSG:26910")
task3.start()