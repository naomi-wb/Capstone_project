## This script is used to loop through a list of species and calculate
## the summary data by month from eBird to then correlate with abundances
## it calls the function from "R/function_to_summarize_eBird_data_for_Birdlife.R"

library(readr)
library(dplyr)

# source function
source("R/Extracting training data from eBird/function_to_summarize_eBird_data_for_Birdlife.R")

# get potential list of species from birdlife file
list_of_species_with_data <- pop_estimates %>%
  dplyr::filter(complete.cases(Lower)) %>%
  dplyr::filter(Lower != "X") %>%
  .$COMMON_NAME

# get list of species which actually worked
list_of_species_with_extracted_data <- data.frame(files=list.files("eBird data/species_dat")) %>%
  mutate(files=as.character(files)) %>%
  mutate(Birdlife=endsWith(.$files, "Birdlife_dat.RDS")) %>%
  dplyr::filter(Birdlife=="TRUE") %>%
  mutate(species_name=gsub("_Birdlife_dat.RDS", "", .$files)) %>%
  mutate(species_name=gsub("_", " ", .$species_name)) %>%
  .$species_name

# Now get nrow (length of data) for each of these datasets
# A lot of them are very remote with little to no data
# so we will not be able to do anything with them
library(magicfor)   

magic_for(put, silent = TRUE)

for (i in unique(list_of_species_with_extracted_data)){
  
  dat <- readRDS(paste0("eBird data/species_dat/", gsub(" ", "_", i), "_Birdlife_dat.RDS"))
  
  number_of_total_obs <- as.numeric(nrow(dat))
  number_of_species_obs <- as.numeric(dat %>%
    dplyr::filter(COMMON_NAME == i) %>%
    nrow(.))
  number_of_checklists <- as.numeric(length(unique(dat$SAMPLING_EVENT_IDENTIFIER)))
                                      
  put(number_of_total_obs,
      number_of_species_obs,
      number_of_checklists)
}

extracted_data_summary <- magic_result_as_dataframe() %>%
  rename(species=i)

# now lets run the function only for species ranges who
# some species give errors, so I manually extracted theme here in the list
# have at least 50 checklists
species_to_analyze <- extracted_data_summary %>%
  dplyr::filter(number_of_checklists > 50) %>%
  dplyr::filter(species != "Cebu Boobook") %>%
  dplyr::filter(species != "Nava's Wren") %>%
  .$species

# for loop to loop through potential species
for (i in unique(species_to_analyze)) {
  
  
  summarize_abundances_for_a_species(i)
  
  
}
