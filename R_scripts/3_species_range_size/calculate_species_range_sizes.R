## Calculate species' range sizes based on
## Their split geojsons
## this will help limit the number of species which should be checked for population
## sizes
################### THIS SCRIPT WILL NOT RUN AS THE SPECIES GEOJSONS ARE NOT AVAILABLE IN THIS REPOSITORY
################### THEY ARE AVAILABLE FROM BirdLife International if you use the script "split_species_range_maps_to_geojsons.R"

# packages
library(sf)
library(dplyr)
library(readr)

# read in clements data
clements <- read_csv("Data/Clements-Checklist-v2018-August-2018.csv") %>%
  dplyr::filter(category == "species") %>%
  rename(COMMON_NAME = `English name`) %>%
  rename(SCIENTIFIC_NAME = `scientific name`) %>%
  dplyr::select(COMMON_NAME, SCIENTIFIC_NAME, order, family)


file_list <- list.files("Data/Species range maps/species_geojsons")


measure_area_function <- function(species_file) {

dat <- st_read(paste0("Data/Species range maps/species_geojsons/", species_file))

area <- as.numeric(st_area(dat))/2.59e+6

sci_name <- gsub(".geojson", "", species_file)
sci_name <- gsub("_", " ", sci_name)

common_name <- clements %>%
  dplyr::filter(SCIENTIFIC_NAME == sci_name) %>%
  .$COMMON_NAME

df <- data.frame(COMMON_NAME=common_name,
                 area_square_miles=area)

return(df)

}

list_of_results <- lapply(file_list, function(x){measure_area_function(x)})

df_results <- bind_rows(list_of_results)

write_csv(df_results, "Data/Species range maps/range_sizes.csv")

# subset to any species with range sizes < 60k square miles
# which is roughly the median of our current states areas values (NA and GB)
# need to minimize the amount of data extracted within a species range
# and minimize the number of species we need to manually search through
potential_species <- df_results %>%
  dplyr::filter(area_square_miles < 60000)

write_csv(potential_species, "Data/Species range maps/potential_birdlife_species_for_inclusion.csv")

