##   ----------------- Extracting Environmental Variables ----------------------

# This script is designed to import modern environmental variables for a given list of  countries.  This results in a metadata file be joined with the sample names (see 7_collect_sequence_names.Rmd) to create a pair of files that can be used to identify correlations and run GLMM analysis. 

#    ---------------------------------------------------------------------------

##   Libraries required ----

if (!require("geodata")) install.packages("geodata")
if (!require("tidyterra")) install.packages("tidyterra")
library(terra)
library(geodata)
library(ggplot2)
library(tidyterra)
library(feather)

## Load in sequence names 
seq_names <- feather::read_feather(
  path = "data/vcf/SequenceNames.txt"
)

##  Importing climate data 

# This downloads data from the world Clim data page 
# [https://www.worldclim.org/data/bioclim.html]

global_bio <- geodata::worldclim_global(
  var = "bio",
  res = 2.5,  # 2.5 minutes resolution
  path = tempdir(), 
  version = "2.1"
)
global_bio
# bio1 : annual mean temp
# bio2 : mean diurnal range
# bio3 : Isothermality (bio2 * bio7)
# bio4 : temperature seasonality
# bio5 : Max temp warmest month 
# bio6 : Min temp coldest month
# bio7 : Temperature annual range 
# bio12: Annual precipitation
# bio15: Precipitation seasonality

# To plot one of these variables...
plot(global_bio[[3]])

# list bioclimatic variables
env_covariates <- global_bio[[c(1,2,3,4,5,6,7,12,15)]]
names(env_covariates) <- c(
  "tmp_yr", "tmp_range_drl", "isotherm", "tmp_seasonality", "tmp_max", "tmp_min", 
  "tmp_range_yr", "precip_yr", "precip_seasonality"
)

# get list of countries 
country_codes <- geodata::country_codes(
  query = NULL
)
# replace name that has changed recently
country_codes$NAME <- gsub(
  pattern = "Czech Republic", 
  replacement = "Czechia", 
  x = country_codes$NAME
)

## generating world coordinates
world_coord <- geodata::world(
  resolution = 3, 
  level = 0, 
  path = tempdir()
)

# extracting averages per country region ----------
env <- terra::extract(
  env_covariates, 
  world_coord, 
  mean,
  na.rm = TRUE, 
  weights = TRUE
)

# convert to dataframe
world_coord_df <- as.data.frame(world_coord$NAME_0)

## combine the mean values with country names...
env <- base::cbind(
  world_coord_df, 
  env
)

# collect country centriod coordinates
world_centroid <- terra::geom(country_centroids)

# combine with mean value data 
env <- base::cbind(
  env, 
  world_centroid
)

env <- env[,c(1,2,15,16,4:12)]
names(env)[1:4] <- c("ISO_3", "country", "lat", "long")


# using names from above...
meta2 <- dplyr::left_join(
  x = seq_names, 
  y = env,
  by = dplyr::join_by(country)
)

# add in regions and continents: 
meta2 <- dplyr::left_join(
  x = meta2, 
  y = country_codes, # quereied above
  by = dplyr::join_by("ISO_3" == "ISO3")
)

# reorder and remove unwanted cols: 
meta2 <- meta2[,c(1:3,16,21,23,4:14)]

# unite first two columns to reform sample names: 
meta2 <- tidyr::unite(
  data = meta2, 
  col = "samples", 
  acc, country,
  sep = "_", 
  na.rm = TRUE, 
  remove = FALSE
)

# remove spaces from sample names: 
meta2$samples <- gsub(
  pattern = " ", 
  replacement = "_", 
  x = meta2$samples
)

# change iso column names (to remove capital letters)
names(meta2)[4:6] <- c("iso3", "iso2", "subcontinent")

# write to alternate file 
feather::write_feather(
  x = meta2, 
  path = "data/meta/worldclim_bioclimatic_data_09_2023.feather"
)

