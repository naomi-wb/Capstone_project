## This script reads in a large geodatabase provided by BirdLife International
## and does some subsampling to then write out each file separately as a geojson
## Then I'll hopefully be able to use the geojson to interact with the ebird bigquery database
## and pull out all records within a geojson
################### THIS SCRIPT WILL NOT RUN AS THE SPECIES range maps (stored as a geodatabase)
################### THEY ARE AVAILABLE FROM BirdLife International

# packages
library(sf)
library(dplyr)
library(readr)

# read in large geodatabase -- takes a few minutes
shape <- st_read("Data/Species range maps/BOTW/BOTW.gdb")

# read in clements data as well
clements <- read_csv("Data/Clements-Checklist-v2018-August-2018.csv") %>%
  dplyr::filter(category == "species") %>%
  rename(COMMON_NAME = `English name`) %>%
  rename(SCIENTIFIC_NAME = `scientific name`) %>%
  dplyr::select(COMMON_NAME, SCIENTIFIC_NAME, order, family)

# filter potential species
shape2 <- shape %>%
  # only resident bird species
  dplyr::filter(SEASONAL == 1) %>%
  # only extant species
  dplyr::filter(PRESENCE == 1) %>%
  # only native species
  dplyr::filter(ORIGIN == 1)

# drop geometry for some summarizing
shape3 <- shape2
st_geometry(shape3) <- NULL

# list of species to make a geojson for
species_list <- shape3 %>%
  group_by(SCINAME) %>%
  summarise(N=n()) %>%
  rename(SCIENTIFIC_NAME = SCINAME) %>%
  inner_join(., clements, by="SCIENTIFIC_NAME") %>%
  dplyr::filter(N==1)

for (i in unique(species_list$SCIENTIFIC_NAME)) {
  
  temp_dat <- shape2 %>%
    dplyr::filter(SCINAME == i)
  
  st_write(temp_dat, paste0("Data/Species range maps/species_geojsons/", gsub(" ", "_", i), ".geojson"))
}









