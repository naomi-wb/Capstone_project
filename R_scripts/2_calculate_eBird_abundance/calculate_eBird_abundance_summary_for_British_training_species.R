# This script is used to loop through a list of species and calculate
## the summary data by month from eBird to then correlate with abundances
## it calls the function from "R/function_to_summarize_eBird_data_for_GB.R"

library(readr)
library(dplyr)

# source function
source("R/Extracting training data from eBird/function_to_summarize_eBird_data_for_GB.R")


## read in cleaned pop estimates
pop_estimates <- read_csv("Data/British Birds data/british_population_estimates.csv")

effort_dat1 <- readRDS("eBird data/effort_dat/GB-NIR.RDS")
effort_dat2 <- readRDS("eBird data/effort_dat/GB-WLS.RDS")
effort_dat3 <- readRDS("eBird data/effort_dat/GB-SCT.RDS")
effort_dat4 <- readRDS("eBird data/effort_dat/GB-ENG.RDS")

effort_dat <- bind_rows(effort_dat1, effort_dat2,
                        effort_dat3, effort_dat4)

rm(effort_dat1)
rm(effort_dat2)
rm(effort_dat3)
rm(effort_dat4)

# for loop to loop through potential species
for (i in unique(pop_estimates$COMMON_NAME)) {
  
  
  summarize_abundances_for_a_species(i)
  
  
}
