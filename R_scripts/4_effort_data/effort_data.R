### An R script to extract the
### effort data
### for every state
### i.e., every unique eBird checklist information
### and pull from bigquery and chunk it into RDS files
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


chunk_state_effort <- function(state) {
  
  state_effort <- ebird %>%
    dplyr::filter(STATE_CODE == state) %>%
    dplyr::select(SAMPLING_EVENT_IDENTIFIER, STATE_CODE, COUNTY, BCR_CODE,
                  LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE,
                  TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                  EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER) %>%
    distinct() %>%
    dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
    collect(n=Inf)
  
  saveRDS(state_effort, paste0("eBird data/effort_dat/", state, ".RDS"))
  
  
}

# first do it for the United States, so will rely on an external
# file which has all the states and the areas of each
# read in list of state codes for US and canada
state_codes <- read_csv("Data/states.csv")

list_of_states <- state_codes %>%
  dplyr::filter(!STATE_CODE %in% c("GB-WLS", "GB-ENG", "GB-NIR", "GB-SCT")) %>%
  .$STATE_CODE


# apply function for every state
# in North America (i.e., USA and Canada)
lapply(list_of_states, function(x) {chunk_state_effort(x)})

# apply function for every state
# in British estimates
list_of_states2 <- state_codes %>%
  dplyr::filter(STATE_CODE %in% c("GB-WLS", "GB-ENG", "GB-NIR", "GB-SCT")) %>%
  .$STATE_CODE

lapply(list_of_states2, function(x) {chunk_state_effort(x)})




