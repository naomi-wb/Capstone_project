## This script is a function to extract data from 
## big query and save it as an RDS for each species
## first need to get a list of species
## from the PIF estimates
## as this is the list of species I'll need to do this for first
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

# get a list of species
# to get data for
PIF_estimates <- read_csv("Data/PIF estimates/PopEsts_BCR_x_ProvState_2-27-19.csv")

colnames(PIF_estimates) <- gsub(" ", "_", colnames(PIF_estimates))

cleaned_pop_estimates <- PIF_estimates %>%
  dplyr::select(1:10) %>%
  replace_na(list(Population_Estimate = "NONE")) %>%
  rename(PROV_STATE = `Province_/_State_/_Territory`) %>%
  dplyr::filter(Population_Estimate != "NONE") %>%
  mutate(Population_Estimate = as.numeric(as.character(Population_Estimate))) 

write_csv(cleaned_pop_estimates, "Data/PIF estimates/cleaned_pop_estimates.csv")


list_of_species <- cleaned_pop_estimates %>%
  distinct(English_Name, .keep_all = TRUE) %>%
  dplyr::select(English_Name, Scientific_Name)


write_csv(list_of_species, "Data/PIF estimates/list_of_species.csv")


# first see which of the list of species exist in the eBird dataset
# with a perfect name match
# will then need to fix the ones which don't match
dat <- ebird %>%
  dplyr::select(COMMON_NAME) %>%
  dplyr::filter(COMMON_NAME %in% list_of_species$English_Name) %>%
  group_by(COMMON_NAME) %>%
  summarise(N=n()) %>%
  collect(n=Inf)

# it looks like three species
# are 'missing' from eBird but are contained in the PIF dataset
# will need to fix these species
list_of_species_corrected <- list_of_species %>%
  rename(COMMON_NAME = English_Name) %>%
  left_join(., dat, by="COMMON_NAME")

# upon further inspection, Northern Hawk Owl, Great Gray Owl, and Gryfalcon are all missing
# these are 'sensitive species' and thus are expected to be missing
# so will ignore this for now and use the 'list of species' as the actual list of species.

# write a function that for each species saves all the data
# as an RDS in a folder
chunk_species_data <- function(species) {
  
  species_dat <- ebird %>%
    dplyr::filter(COMMON_NAME == species) %>%
    dplyr::select(SAMPLING_EVENT_IDENTIFIER, STATE_CODE, COUNTY, BCR_CODE,
                  LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE,
                  TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                  EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER, OBSERVATION_COUNT) %>%
    dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
    collect(n=Inf)
  
  species2 <- gsub(" ", "_", species)
  
  saveRDS(species_dat, paste0("eBird data/species_dat/", species2, ".RDS"))
  
  
}


lapply(list_of_species %>%
         dplyr::filter(English_Name != "Northern Hawk Owl") %>%
         dplyr::filter(English_Name != "Gyrfalcon") %>%
         dplyr::filter(English_Name != "Great Gray Owl") %>%
         .$English_Name,
       function(x) {chunk_species_data(x)})


