# this script is a script to
# extract data for birdlife population estimates to include in the analysis
################ THIS SCRIPT WILL NOT WORK AS IT REQUIRES ACCESS TO GOOGLE BIG QUERY WHERE THE EBIRD DATA ARE STORED


## packages
library(readr)
library(bigrquery)
library(dbplyr)
library(dplyr)
library(tidyr)
library(sf)
library(geojsonsf)
library(DBI)


# read in clements data
clements <- read_csv("Data/Clements-Checklist-v2018-August-2018.csv") %>%
  dplyr::filter(category == "species") %>%
  rename(COMMON_NAME = `English name`) %>%
  rename(SCIENTIFIC_NAME = `scientific name`) %>%
  dplyr::select(COMMON_NAME, SCIENTIFIC_NAME, order, family)

# read in assessed birdlife species
# and join with clements so that scientific name is also available
birdlife_species <- read_csv("Data/Species range maps/ASSESSED_potential_birdlife_species_for_inclusion.csv") %>%
  left_join(., clements, by="COMMON_NAME")

# create connection with online database
con <- DBI::dbConnect(bigrquery::bigquery(),
                      dataset= "ebird",
                      project="ebird-database",
                      billing="ebird-database")

# function to extract dat for a species in a list
extract_species_dat_for_Birdlife <- function(species){

# filter the data for the species
species_sliced <- birdlife_species %>%
  dplyr::filter(COMMON_NAME == species)

# read in geojson range
range <- st_read(paste0("Data/Species range maps/species_geojsons/", gsub(" ", "_", species_sliced$SCIENTIFIC_NAME), ".geojson"))

# get geojson string from the sf object
geojson_string <- sf_geojson(range)

# manipulate the string to put into
# a paste0 below to extract from BigQuery
string <- gsub('.{3}$', '', sub('.*coordinates', '', geojson_string))

# extract all dat from within the geojson string
species_dat <- dbGetQuery(con, 
                    paste0("SELECT
                           SAMPLING_EVENT_IDENTIFIER, STATE_CODE, LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE,
                           TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, COMMON_NAME,
                           SCIENTIFIC_NAME, EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, CATEGORY,
                           GROUP_IDENTIFIER, OBSERVATION_COUNT
                           FROM
                           ebird.ebird_qa
                           WHERE
                           ST_CONTAINS(ST_GEOGFROMGEOJSON('{\"type\":\"MultiPolygon\",\"coordinates", string,
                           "'),
                           ST_GEOGPOINT(LONGITUDE, LATITUDE))")
                    )

saveRDS(species_dat, paste0("eBird data/species_dat/", gsub(" ", "_", species), "_Birdlife_dat.RDS"))

}

list_of_species_with_data <- birdlife_species %>%
  dplyr::filter(complete.cases(Estimate)) %>%
  dplyr::filter(Lower != "X") %>%
  .$COMMON_NAME


lapply(list_of_species_with_data, function(x) {extract_species_dat_for_Birdlife(x)})

FaultTolerantFunction <- function(species) {
tryCatch({ret <- extract_species_dat_for_Birdlife(species);}, error = function(e) {ret <<- NA});
ret
}

lapply(list_of_species_with_data, function(x) {FaultTolerantFunction(x)})

