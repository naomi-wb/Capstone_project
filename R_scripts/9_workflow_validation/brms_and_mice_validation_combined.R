# brms and mice validation
# loading packages
pacman::p_load(tidyverse, 
               purrr,
               magrittr,
               MCMCglmm,
               broom,
               lme4,
               ggplot2,
               patchwork,
               scales,
               rstan,
               brms,
               tidybayes,
               forcats,
               readr,
               GGally,
               corrplot,               
               here,
               rapportools,
               mvtnorm,
               sf,
               ggcorrplot,
               mice, 
               miceadds, 
               mitml
)


# running rstan + brms
options(mc.cores = 24)
rstan_options(auto_write = TRUE)

# read in the data for detectability modelling
mod_df <- readRDS("Intermediate data files/data_for_detectability_modelling.RDS")

# read in data that was imputed
# necessary for the factors and such within
dat <- readRDS("Intermediate data files/file_for_imputation.RDS")

# read in other data that will be used later
# read in clements clean
clements_clean <- read_csv("Data/clements_clean.csv")

# read in body size clean
body_size <- read_csv("Data/body size data/cleaned_body_size_data.csv") %>%
  dplyr::select(3, 7:9)

# read in color data
color <- read_csv("Data/color data/color_data_by_species.csv")

# IUCN data
IUCN <- read_csv("Data/IUCN categories/cleaned_IUCN_data.csv")

# flock size
flock_size <- readRDS("Data/Flock size/flock_size_per_month_for_each_species.RDS") %>%
  group_by(COMMON_NAME) %>%
  summarise(mean_max_flock_size=mean(max_abund)) %>%
  dplyr::rename(ebird_COMMON_NAME=COMMON_NAME)

# read in the data files, which were stratified by month
# to make storage of them easier
# and in case we want to stratify analyses to different months
january <- readRDS("Intermediate data files/abundance_grid_results/january.RDS")
february <- readRDS("Intermediate data files/abundance_grid_results/february.RDS")
march <- readRDS("Intermediate data files/abundance_grid_results/march.RDS")
april <- readRDS("Intermediate data files/abundance_grid_results/april.RDS")
may <- readRDS("Intermediate data files/abundance_grid_results/may.RDS")
june <- readRDS("Intermediate data files/abundance_grid_results/june.RDS")
july <- readRDS("Intermediate data files/abundance_grid_results/july.RDS")
august <- readRDS("Intermediate data files/abundance_grid_results/august.RDS")
september <- readRDS("Intermediate data files/abundance_grid_results/september.RDS")
october <- readRDS("Intermediate data files/abundance_grid_results/october.RDS")
november <- readRDS("Intermediate data files/abundance_grid_results/november.RDS")
december <- readRDS("Intermediate data files/abundance_grid_results/december.RDS")

# combine these data
combined_dat <- bind_rows(january, february, march,
                          april, may, june, july, 
                          august, september, october,
                          november, december)

# read in grid cell land mass estimates
grids_land_area <- readRDS("Data/Spatial_data/grid_percent_land.RDS") %>%
  dplyr::select(ID, percentage) %>%
  dplyr::rename(grid_id=ID) %>%
  st_set_geometry(NULL)

area_summaries <- readRDS("Intermediate data files/additional_area_summary.RDS")

# get the total area (by grids) that a species occupies
# and join with the range maps area
# and adjust the areas accordingly
total_area_per_species <- dat %>%
  group_by(ebird_COMMON_NAME) %>%
  summarize(total_area_grids=sum(area_square_miles)) %>%
  left_join(., area_summaries) %>%
  mutate(total_area=ifelse(range_area_square_miles < 26000 & range_area_square_miles > 0, 
                           range_area_square_miles, total_area_grids)) %>%
  mutate(total_area=ifelse(diff_area_square_miles>0, total_area_grids+diff_area_square_miles, total_area)) %>%
  mutate(total_area=ifelse(is.na(total_area)==TRUE, total_area_grids, total_area)) %>%
  mutate(range_adjusted=ifelse(total_area_grids==total_area, FALSE, TRUE))

# checking for infinity
which(mod_df$density_log10 == -Inf)
m_inf <- which(mod_df$abundance_log10 == -Inf)
min(mod_df$abundance_log10[-m_inf]) # the smallest is -4.499787
# checking for NA 
which(is.na(mod_df$abundance_SD_log10) == T) # SD = SE is not avaiable
max(mod_df$abundance_SD_log10, na.rm = T) # the largest is 0.5063101

# so we replace -Inf with -4.5 for abundance_log10
# there is one count == 1 for count_sd - making it 2 so we can calcuate SE for aboundance_log10
# we use sample size corrected 
mod_df2 <- mutate(mod_df, abundance_log10 = if_else(abundance_log10 == -Inf, -4.5, abundance_log10),
                  abundance_SD_log10 = if_else(is.na(abundance_SD_log10) == T, 0.51, abundance_SD_log10)) 

# Now a function to take a species
# withold that species
# run the brms without that species
# run the imputation without that species
# and then summarize the global population for that species
# THEN see if when that species was a training species
# and when it wasn't a training species are correlated
# don't need to do this for too many species I don't think.

# re-do analysis function
re_do_analysis <- function(species_name){
  
  message(paste0("Analyzing", species_name))
  
  mod_dat <- mod_df2 %>%
    dplyr::filter(COMMON_NAME!=species_name)
  
  
  brms_mod_final <- brm(density_log10 ~ 1 + me(abundance_log10, abundance_SD_log10) + 
                          (1 + me(abundance_log10, abundance_SD_log10) |species), 
                        data = mod_dat, save_mevars = T, chains = 2, iter = 3000, warmup = 1000)
  
  
  # Summarize data for this model
  means <- fixef(brms_mod_final)
  blups<-ranef(brms_mod_final)
  
  intercepts <- means[1,1] + as.numeric(blups$species[ , 1, 1])
  slopes <-   means[2,1] + + as.numeric(blups$species[ , 1, 2])
  
  # errors
  intercepts_se <- sqrt(means[1,2]^2 + as.numeric(blups$species[ , 2, 1])^2)
  slopes_se <-   sqrt(means[2,2]^2 + as.numeric(blups$species[ , 2, 2])^2)
  
  n_spp <- length(unique(mod_dat$species))
  blups_df_final <- tibble(species = rep(attr(blups$species, "dimnames")[[1]], 2),
                           type = rep(c("intercept","slope"), each = n_spp),
                           estimate = c(intercepts, slopes), 
                           se =  c(intercepts_se, slopes_se))
  
  # read in results of MCMCglmm model
  detectability <- blups_df_final %>%
    pivot_wider(names_from=type, values_from=c(estimate, se)) %>%
    dplyr::rename(ebird_COMMON_NAME=species) %>%
    dplyr::rename(intercept=estimate_intercept) %>%
    dplyr::rename(slope=estimate_slope) %>%
    dplyr::rename(intercept_se=se_intercept) %>%
    dplyr::rename(slope_se=se_slope)
  
  # Now read in all the data used for imputation
  # and create a clean file of these data
  imputation_file <- clements_clean %>%
    distinct(ebird_COMMON_NAME, .keep_all=TRUE) %>%
    dplyr::select(ebird_COMMON_NAME, TipLabel) %>%
    left_join(., flock_size, by="ebird_COMMON_NAME") %>%
    left_join(., color, by="TipLabel") %>%
    left_join(., body_size, by="TipLabel") %>%
    left_join(., IUCN, by="TipLabel") %>%
    group_by(ebird_COMMON_NAME) %>%
    slice(1) %>%
    # clean the data based on IUCN categories
    # remove any species which are 'Extinct in the wild"
    mutate(IUCN_ordinal = case_when(
      IUCN_category == "EW" ~ 0,
      IUCN_category == "EX" ~ 0,
      IUCN_category == "CR (PE)" ~ 0,
      IUCN_category == "CR (PEW)" ~ 0,
      IUCN_category == "CR" ~ 1,
      IUCN_category == "EN" ~ 2,
      IUCN_category == "VU" ~ 3,
      IUCN_category == "NT" ~ 4,
      IUCN_category == "LC" ~ 5
    )) %>%
    left_join(., detectability) %>%
    dplyr::select(1:8, 13:18) %>%
    ungroup() %>%
    mutate(adult_body_mass_g=ifelse(adult_body_mass_g==-999, NA, adult_body_mass_g)) %>%
    left_join(., clements_clean %>%
                dplyr::select(ebird_COMMON_NAME, ebird_SCIENTIFIC_NAME,
                              order, family) %>%
                distinct())
  
  temporal_summary <- combined_dat %>%
    dplyr::filter(mean_abund>0) %>%
    group_by(grid_id, ebird_COMMON_NAME) %>%
    summarize(mean_abund=mean(mean_abund),
              number_months=length(unique(MONTH)),
              mean_species_checklists=mean(obs_with_abund),
              total_species_checklists=sum(obs_with_abund)) %>%
    ungroup() %>%
    left_join(., combined_dat %>%
                dplyr::filter(mean_abund>0) %>%
                dplyr::select(grid_id, total_checklists) %>%
                distinct() %>%
                group_by(grid_id) %>%
                summarize(mean_grid_checklists=mean(total_checklists),
                          total_grid_checklists=sum(total_checklists)), by="grid_id") %>%
    left_join(., combined_dat %>%
                dplyr::select(grid_id, area_square_miles) %>%
                distinct(), by="grid_id")
  
  # only select species with a positive mean abundance
  # assumes that species with 0 mean abundance are 'true' zeros
  dat <- temporal_summary 
  
  # transforming data
  imputation_file %<>% mutate(mass_log10 = log10(adult_body_mass_g),
                              flock_size_log10 = log10(mean_max_flock_size),
                              species = ebird_COMMON_NAME) 
  
  dat_imp <- imputation_file %>% 
    dplyr::select(ebird_COMMON_NAME, max_distance, max_brightness,
                  mass_log10, flock_size_log10, IUCN_ordinal) %>% 
    as.data.frame()
  
  # join things together
  joined_dat <- dat %>%
    left_join(., dat_imp) %>%
    left_join(., detectability) %>%
    mutate(density=10^(intercept+slope*(log10(mean_abund)))) %>%
    mutate(log10_density=log10(density)) %>%
    mutate(log10_mean_abund=log10(mean_abund))
  
  # Now need to calculate the SE for every density estimate that we have
  getSE <- function(species_name, cor =  -1, sim_N = 10000){
    
    spp <- joined_dat %>%
      dplyr::filter(ebird_COMMON_NAME == species_name)
    
    # now add a subsection to apply this across every grid
    # that a species occurs in
    grid_split_function <- function(id_of_grid, species_name) {
      
      dat <- spp %>%
        dplyr::filter(grid_id == id_of_grid)
      
      # prep for our varaince-covariance matrix
      cov <- dat$slope_se*dat$intercept_se*cor
      int_v <- dat$intercept_se^2  
      slo_v <- dat$slope_se^2
      
      vcv <- matrix(c(int_v, cov, cov, slo_v), nr = 2)
      beta <- c(dat$intercept,dat$slope)
      
      # 10000s of the intercept-slope pairs
      beta_sim <- rmvnorm(sim_N, beta, vcv)
      
      # creating the desgin matrix
      x <-cbind(1, dat$log10_mean_abund)
      
      # some matrix algebra
      xb <- x %*% t(beta_sim)
      
      predicted_density <- apply(xb,1,quantile, probs = c(.5,0.8413448)) %>%
        t() %>%
        as.data.frame()
      
      se <- data.frame(grid_id=id_of_grid,
                       se=predicted_density[,2]-predicted_density[,1])
      
      return(se)
    }
    
    # apply over grids
    species_summary <- bind_rows(lapply(unique(spp$grid_id), grid_split_function)) %>%
      mutate(ebird_COMMON_NAME=species_name)
    
    return(species_summary)
    
  }
  
  # split joined dat
  # to those that have density estimates (need to get SE for)
  # and those who don't
  training_dat <- joined_dat %>%
    dplyr::filter(complete.cases(log10_density))
  
  length(unique(training_dat$ebird_COMMON_NAME))
  
  data_for_imputation <- joined_dat %>%
    dplyr::filter(is.na(log10_density))
  
  # apply the function to calculate SE for every species/row in the training data
  SEs_1 <- bind_rows(lapply(unique(training_dat$ebird_COMMON_NAME), function(x){getSE(x, cor=-1)})) %>%
    mutate(Cor=-1)
  
  training_dat.2 <- training_dat %>%
    left_join(., SEs_1) %>%
    dplyr::select(-Cor)
  
  data_for_imputation.2 <- data_for_imputation %>%
    mutate(se=NA)
  
  joined_dat_final <- bind_rows(training_dat.2, data_for_imputation.2) %>%
    left_join(., grids_land_area)
  
  
  sdat <- joined_dat_final %>% 
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
              m = 10, # need to get at least 5 - perferably 10
              maxit = 20, # we probably neeed to 20 to coverge
              method = imp_method,
              predictorMatrix = pred_matrix,
              seed = 777)
  
  # get grid and species level checklists
  # as factors to join with imputation data
  factors <- joined_dat_final %>%
    dplyr::select(ebird_COMMON_NAME, grid_id, 
                  total_grid_checklists, total_species_checklists) %>%
    ungroup() %>%
    mutate(`.id`=1:nrow(.))
  
  # get grid and species level checklists
  grid_species_checklists <- joined_dat_final %>%
    dplyr::select(ebird_COMMON_NAME, grid_id, 
                  total_grid_checklists, total_species_checklists)
  
  # get imputed data as a list
  imp_list <- mice::complete(imp, "long")
  
  # join factors with imputation results
  idat <- imp_list %>% 
    left_join(., factors)
  
  #  collapsing imputations
  # so taking values across the 100 imputed datasets
  idat2 <- idat %>% 
    group_by(ebird_COMMON_NAME, grid_id) %>% 
    summarise(imp_se = sd(density, na.rm = T), 
              density = mean(density),
              dens_se = mean(dens_se), 
              se = sqrt(imp_se^2 + dens_se^2),
              se = ifelse(se > 1, 1, se)) # I am making the max of SE = 1
  
  # summarize with confidence intervals
  ddat <- idat2 %>% 
    mutate(mode = exp((density*log(10) - (se*log(10))^2)),
           median = exp(density*log(10)),
           mean = exp((density*log(10) + 0.5*(se*log(10))^2)),
           CILL = exp((density - se*qnorm(0.975))*log(10)), # 95% CI on log10 scale turn onto original
           CIUL = exp((density + se*qnorm(0.975))*log(10))) # this already include imputation SE (which now we trancate at 1)
  
  # investigate the relationship between density
  # and the number of checklists for a species in a grid
  # and the total number of grid checklists
  checklist_dens <- ddat %>% 
    ungroup () %>% 
    mutate(grid_id = as.numeric(as.character(grid_id))) %>%
    group_by(ebird_COMMON_NAME) %>%
    left_join(., grid_species_checklists) %>%
    ungroup()
  
  # weighted se function
  weighted.ave.se <- function(se, weights) { # se = vector, weights = the number of lists
    if(length(se) != length(weights)) stop('se and weights are not the same length')
    var <- se^2
    k <- length(se)
    w.var<-sum(var*(weights-1))/(sum(weights)-k)
    w.se <- sqrt(w.var)
    return(w.se)
  }
  
  # weight the densities by
  # the number of checklists
  # providing more confidence to species with more checklists
  weighted_dens <- ddat %>% 
    ungroup () %>% 
    mutate(grid_id = as.numeric(as.character(grid_id))) %>%
    group_by(ebird_COMMON_NAME) %>%
    left_join(., grid_species_checklists) %>%
    mutate(weights=total_species_checklists/total_grid_checklists) %>%
    #mutate(weights=1) %>%
    summarise(density = weighted.mean(density, weights),
              se_non = sqrt(weighted.ave.se(dens_se^2, weights)),
              se_imp = sqrt(weighted.ave.se(se^2, weights)))
  
  # now get abundance estimates
  sppdat <- weighted_dens %>%
    left_join(., total_area_per_species)
  
  # abundance for all the species as ditribuitons of 10,000
  list_dens_non <- map2(sppdat$density, sppdat$se_non, ~ rlnorm(1e4, .x*log(10), .y*log(10))) %>% 
    purrr::map(~.[order(.)]) 
  list_dens_imp <- map2(sppdat$density, sppdat$se_imp, ~ rlnorm(1e4, .x*log(10), .y*log(10))) %>% 
    purrr::map(~.[order(.)])
  
  # 10,000 draws for all the spp 
  spp_num_non <- map2_dfc(list_dens_non, sppdat$total_area, ~.x*.y)              
  spp_num_imp <- map2_dfc(list_dens_imp, sppdat$total_area, ~.x*.y)  
  
  # summing these draws
  # summing is the step - means between two things become different
  spp_sum_non <- as.vector(rowSums(spp_num_non))
  spp_sum_imp <- as.vector(rowSums(spp_num_imp))
  
  # quantile = confidence intervals; hdi = higest desnity (much like Bayesian higest posterior density)
  species_col_number <- sppdat %>%
    mutate(row_id=1:nrow(.)) %>%
    dplyr::filter(ebird_COMMON_NAME == species_name) %>%
    .$row_id
  
  species_dat <- tibble(non = spp_num_non[[species_col_number]],
                        imp = spp_num_imp[[species_col_number]]) %>%
    mutate(ebird_COMMON_NAME=sppdat$ebird_COMMON_NAME[[species_col_number]])
  
  saveRDS(species_dat, paste0("Intermediate data files/workflow_validation/with_held/", species_name, ".RDS"))
}


# get a list of the observed species
# to repeat the imputation function
# above for
# randomly sample this for 100 species now
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

obs_species <- factors %>%
  dplyr::filter(Data=="observed") %>%
  dplyr::select(ebird_COMMON_NAME) %>%
  distinct() %>%
  .$ebird_COMMON_NAME

mclapply(obs_species, re_do_analysis)




