## This script gets the data from eBird for british birds
## from 'Great Britain' state codes
################ THIS SCRIPT WILL NOT WORK AS IT REQUIRES ACCESS TO GOOGLE BIG QUERY WHERE THE EBIRD DATA ARE STORED

## packages
library(readr)
library(bigrquery)
library(dbplyr)
library(dplyr)
library(tidyr)


# create connection with online database
con <- DBI::dbConnect(bigrquery::bigquery(),
                      dataset= "ebird",
                      project="ebird-database",
                      billing="ebird-database")

# create ebird table
ebird <- tbl(con, 'ebird_qa')

# read in british data
british_data <- read_csv("Data/British Birds data/british_population_estimates.csv")

list_of_species <- british_data %>%
  distinct(COMMON_NAME) %>%
  .$COMMON_NAME


chunk_species_data <- function(species) {
  
  species_dat <- ebird %>%
    dplyr::filter(STATE_CODE %in% c("GB-NIR", "GB-SCT", "GB-WLS", "GB-ENG")) %>%
    dplyr::filter(COMMON_NAME == species) %>%
    dplyr::select(SAMPLING_EVENT_IDENTIFIER, STATE_CODE, LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE,
                  TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                  EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER, OBSERVATION_COUNT) %>%
    dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
    collect(n=Inf)
  
  species2 <- gsub(" ", "_", species)
  
  saveRDS(species_dat, paste0("eBird data/species_dat/", species2, "_British_dat.RDS"))
  
}


lapply(list_of_species, function(x) {chunk_species_data(x)})





