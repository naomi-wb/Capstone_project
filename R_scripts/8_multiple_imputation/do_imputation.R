# this is a script to perform the imputations
# at different levels
# and save out the results
# it first gets data from a dataset
# that was developed with the "prepare_file_for_imputation.R" script
# this script runs imputations and saves out the results

# load packages
# gettng packages ready
pacman::p_load(readr, dplyr, lme4, mice, miceadds, 
               mitml, GGally, purrr, HDInterval, scales)


# read in file for imputation
# data with missingness 
dat <- readRDS("Intermediate data files/file_for_imputation_15p.RDS")


# prepare data for imputation
# by logging or standardizing values
sdat <- dat %>% 
  ungroup() %>%  
  dplyr::select(grid_id, total_grid_checklists, area_square_miles,
                ebird_COMMON_NAME, total_species_checklists, number_months, IUCN_ordinal,
                max_distance, max_brightness, mass_log10, flock_size_log10,
                log10_density, log10_mean_abund, se) %>% 
  mutate_if(is.character, as.factor) %>% 
  transmute(species = as.integer(as.factor(ebird_COMMON_NAME)), 
            #grid_id = as.factor(grid_id), # without grid info
            #area = log(area_square_miles),
            #checklistG = log(total_grid_checklists+0.5),
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


# dim to check if the dimention is reduced
dim(sdat)

# look at distribution of se
hist(sdat$dens_se)

md.pattern(sdat)

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

# perform the imputation
imp <- mice(sdat,
            m = 100, # need to get at least 5 - perferably 10
            maxit = 20, # we probably neeed to 20 to coverge
            method = imp_method,
            predictorMatrix = pred_matrix,
            seed = 777)




# write out imputation results
saveRDS(imp, "Intermediate data files/imputation_results_15p.RDS")













