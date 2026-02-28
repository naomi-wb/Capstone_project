## Script is to get a measure of 'flock size'
## for every species
## which is likely to influence the probability of it being detected

## packages
library(readr)
library(bigrquery)
library(dbplyr)
library(dplyr)
library(tidyr)
library(lubridate)


# create connection with online database
con <- DBI::dbConnect(bigrquery::bigquery(),
                      dataset= "ebird",
                      project="ebird-database",
                      billing="ebird-database")

# create ebird table
ebird <- tbl(con, 'ebird_qa')

# testing first
flock_size_per_month <- ebird %>%
  dplyr::select(COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT, 
                CATEGORY, OBSERVATION_DATE) %>%
  dplyr::filter(CATEGORY %in% c("species", "issf", "domestic")) %>%
  dplyr::filter(OBSERVATION_COUNT != "X") %>%
  mutate(OBSERVATION_COUNT=as.numeric(as.character(OBSERVATION_COUNT))) %>%
  mutate(Month = extract(NULL %month from% OBSERVATION_DATE)) %>%
  group_by(COMMON_NAME, SCIENTIFIC_NAME, Month) %>%
  summarise(mean_abund=mean(OBSERVATION_COUNT),
            max_abund=max(OBSERVATION_COUNT),
            median_abund=median(OBSERVATION_COUNT)) %>%
  collect(n=Inf)

saveRDS(flock_size_per_month, "Data/Flock size/flock_size_per_month_for_each_species.RDS")



