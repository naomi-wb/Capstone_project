# mice validation

# packages
library(readr)
library(dplyr)
library(purrr)
library(HDInterval)
library(scales)
library(ggplot2)
library(mice)
library(miceadds)
library(VIM)
library(patchwork)

# read in data that was imputed
# necessary for the factors and such within
dat <- readRDS("Intermediate data files/file_for_imputation.RDS")

# get grid and species level checklists
# as factors to join with imputation data
factors <- dat %>%
  dplyr::select(ebird_COMMON_NAME, grid_id, 
                total_grid_checklists, total_species_checklists,
                intercept) %>%
  ungroup() %>%
  mutate(`.id`=1:nrow(.)) %>%
  mutate(Data=ifelse(is.na(intercept)==TRUE, "imputed", "observed")) %>%
  dplyr::select(-intercept)


# summarize imputation results
imp <- readRDS("Intermediate data files/imputation_results.RDS")

# Now we need to run mice a few times to see how it works
# and for some cross-validation



# A function to apply mice but removing one species
# data points each time
# one of the species from our observed datasets
# And then summarize the information for that species
# that was imputed
# for comparison with the observed values
# the function writes out each species'
# summarized values at a time
mouse_function <- function(species_name){
  
  message(paste0("Analyzing ", species_name))
  
  # prepare data for imputation
  # by logging or standardizing values
  sdat <- dat %>% 
    ungroup() %>%  
    dplyr::select(grid_id, total_grid_checklists, area_square_miles,
                  ebird_COMMON_NAME, total_species_checklists, number_months, IUCN_ordinal,
                  max_distance, max_brightness, mass_log10, flock_size_log10,
                  log10_density, log10_mean_abund, se) %>% 
    mutate_if(is.character, as.factor) %>% 
    mutate(log10_density=ifelse(ebird_COMMON_NAME==species_name, NA, log10_density)) %>%
    mutate(se=ifelse(ebird_COMMON_NAME==species_name, NA, se)) %>%
    transmute(species = as.integer(as.factor(ebird_COMMON_NAME)), 
              checklists = log(total_species_checklists + 0.5),
              months = number_months, 
              IUCN = IUCN_ordinal,
              distance = log10(max_distance),
              brightness = log10(max_brightness),
              mass = mass_log10,
              flock = flock_size_log10,
              abundance = log10_mean_abund,
              density = log10_density,
              dens_se =  se)
  
  
  # set up redictor matrix and imputation methods
  pred_matrix <- make.predictorMatrix(sdat)
  imp_method <- make.method(sdat)
  
  # our cluster (-2)
  pred_matrix[ , "species"] <- -2 # cluster varaible needs to be intergers
  
  # stting 0 for non-missing data
  no_missing <- c("species", "checklists", "months", "abundance")
  pred_matrix[no_missing, ] <- 0
  
  # also put 0 for diag
  diag(pred_matrix) <- 0
  
  pred_matrix
  
  # setting 
  imp_method[c("flock", "IUCN", "mass", "distance", "brightness")] <- "2lonly.pmm" # individual (spp) level 
  
  imp_method[c("density", "dens_se")] <- "2l.pmm" # obs level
  
  imp_method
  
  imp <- mice(sdat,
              m = 10, # running it for low number for model validation stuff only
              maxit = 10, 
              method = imp_method,
              predictorMatrix = pred_matrix,
              seed = 777)
  
  # get imputed data as a list
  imp_list_tmp <- mice::complete(imp, "long")
  
  # join factors with imputation results
  idat_tmp <- imp_list_tmp %>% 
    left_join(., factors)
  
  sp_dat <- idat_tmp %>%
    dplyr::filter(ebird_COMMON_NAME == species_name) %>%
    dplyr::select(1, 11:15) %>%
    mutate(Data="imputation")
  
  saveRDS(sp_dat, paste0("Intermediate data files/mice_validation/", species_name, ".RDS"))
  
}


# get a lit of the observed species
# to repeat the imputation function above
obs_species <- factors %>%
  dplyr::filter(Data=="observed") %>%
  dplyr::select(ebird_COMMON_NAME) %>%
  distinct()

# now apply the function over these species
# will probably take some time!
lapply(obs_species$ebird_COMMON_NAME, mouse_function)


