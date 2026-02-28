## This script is used to loop through a list of species and calculate
## the summary data by month from eBird to then correlate with abundances
## it calls the function from "R/function_to_summarize_eBird_data_for_NA.R"

library(readr)
library(dplyr)

# source function
source("R/Extracting training data from eBird/function_to_summarize_eBird_data_for_NA.R")

## read in cleaned pop estimates for North America
pop_estimates <- read_csv("Data/PIF estimates/cleaned_pop_estimates.csv")


## get a list of species
## but remove three species we don't have data for
list_of_species <- read_csv("Data/PIF estimates/list_of_species.csv") %>%
  dplyr::filter(English_Name != "Northern Hawk Owl") %>%
  dplyr::filter(English_Name != "Gyrfalcon") %>%
  dplyr::filter(English_Name != "Great Gray Owl")





# for loop to loop through potential species
for (i in unique(list_of_species$English_Name)) {
  
  
  summarize_abundances_for_a_species(i)
  
  
}


