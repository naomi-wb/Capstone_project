# This script cleans the extracted abundance data
# and zero-fills the dataframes

## packages
library(readr)
library(dplyr)
library(tidyr)
library(sf)


# will use our clean clements dataset for this
clements <- read_csv("Data/clements_clean.csv") %>%
  group_by(ebird_COMMON_NAME, ebird_SCIENTIFIC_NAME) %>%
  slice(1) %>%
  ungroup()

# detectability data
detectability <- readRDS("Intermediate data files/imputation_results_all.RDS")

# write a function to summarize the data for each
# month
get_grid_abundances_function <- function(month_name){
  
  dat <- readRDS(paste0("Intermediate data files/global_mean_abundance_from_ebird/", month_name, ".RDS"))
  
  y <- rep(1, length(unique(clements$ebird_COMMON_NAME)))
  
  for (i in unique(dat$grid_id)){
    x <- rep(i, 10585)
    y <- c(y, x)
  }
  
  z <- y %>%
    as.data.frame() %>%
    dplyr::rename(grid_id=1) %>%
    dplyr::filter(grid_id != 1) %>%
    .$grid_id
  
  # combine to a df for joining
  full_df_for_joining <- data.frame(COMMON_NAME=rep(unique(clements$ebird_COMMON_NAME), length(unique(dat$grid_id))),
                                    grid_id=z)
  
  dat_joined <- dat %>%
    right_join(., full_df_for_joining) %>%
    ungroup() %>%
    dplyr::select(grid_id, COMMON_NAME, mean_abund, obs_with_abund) %>%
    left_join(., dat %>%
                ungroup() %>%
                dplyr::select(grid_id, area_square_miles, total_checklists) %>%
                distinct(), by="grid_id") %>%
    mutate(MONTH=month_name) %>%
    dplyr::rename(ebird_COMMON_NAME=COMMON_NAME) %>%
    replace_na(list(mean_abund=0))
  
  saveRDS(dat_joined, paste0("Intermediate data files/abundance_grid_results/", month_name, ".RDS"))
}

get_grid_abundances_function("january")
get_grid_abundances_function("february")
get_grid_abundances_function("march")
get_grid_abundances_function("april")
get_grid_abundances_function("may")
get_grid_abundances_function("june")
get_grid_abundances_function("july")
get_grid_abundances_function("august")
get_grid_abundances_function("september")
get_grid_abundances_function("october")
get_grid_abundances_function("november")
get_grid_abundances_function("december")

