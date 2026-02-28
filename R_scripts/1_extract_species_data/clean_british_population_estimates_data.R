## This script is to clean the british birds data more
## I already cleaned it manually in excel
## but this trims it to species we might actually use

# packages
library(readr)
library(dplyr)
library(tidyr)

data <- read_csv("Data/British Birds data/manually_cleaned_british_birds_data.csv")


clean_data <- data %>%
  # only include estimates with a region of Great Britain
  dplyr::filter(Region == "GB") %>%
  # only include estimates with pairs or individuals
  dplyr::filter(Units %in% c("P", "I")) %>%
  # add the midpoint if only a range was given
  mutate(Estimate=case_when(
    is.na(Estimate)=="TRUE" ~ (Lower+Upper)/2,
    is.na(Estimate)=="FALSE" ~ Estimate
  )) %>%
  # only select estimates with reliability 1 or 2
  dplyr::filter(Reliability %in% c(1, 2)) %>%
  # filter complete cases
  dplyr::filter(complete.cases(COMMON_NAME)) %>%
  # convert numbers based on pairs or individuals
  mutate(Estimate=case_when(
    Units == "P" ~ Estimate*2,
    Units == "I" ~ Estimate
  )) %>%
  # repeat the above for lower and upper estimates
  mutate(Lower=case_when(
    Units == "P" ~ Lower*2,
    Units == "I" ~ Lower
  )) %>%
  mutate(Upper=case_when(
    Units == "P" ~ Upper*2,
    Units == "I" ~ Upper
  ))

# export file
write_csv(clean_data, "Data/British birds data/british_population_estimates.csv")
  
  
  
  
  
  
           
          
