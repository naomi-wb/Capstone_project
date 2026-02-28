# This function is used to read in the imputation results
# performed from the script "do_imputation.R"
# and then saved as RDS
# Here we read this in and summarize these and do random draws to 
# derive distribution of possible abundance measures at various levels
# and collapse the imputations accounting for error across each

# packages
library(readr)
library(dplyr)
library(purrr)
library(HDInterval)
library(scales)
library(ggplot2)
library(tidyr)

# read in clean clements
# which will come in handy above
clements <- read_csv("Data/clements_clean.csv")

# read in ecological niche and grouping data
niches <- read_csv("Data/ecological_niche_assignment/41559_2019_1070_MOESM3_ESM.csv") %>%
  dplyr::rename(TipLabel=Binomial)

# read in data that was imputed
# necessary for the factors and such within
dat <- readRDS("Intermediate data files/file_for_imputation_15p.RDS") #using name of new file

# get grid and species level checklists
# as factors to join with imputation data
factors <- dat %>%
  dplyr::select(ebird_COMMON_NAME, grid_id, 
                total_grid_checklists, total_species_checklists) %>%
  ungroup() %>%
  mutate(`.id`=1:nrow(.))

# get grid and species level checklists
grid_species_checklists <- dat %>%
  dplyr::select(ebird_COMMON_NAME, grid_id, 
                total_grid_checklists, total_species_checklists)

# need to get the total area
# for each species where the number needs to be
# multiplied by
# but this needs to be adjusted
# to then get over the Sahara/Siberia problem
# and also drop species' with low area that doesn't encompass a full grid size
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


# summarize imputation results
imp <- readRDS("Intermediate data files/imputation_results_15p.RDS") #adjust new file name

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

# make a plot of the relationship between density and the
# 'weights' in a grid cell used to get a weighted density
ggplot(checklist_dens, aes(y=10^density, x=total_species_checklists/total_grid_checklists))+
  #scale_x_log10(labels=comma)+
  #scale_y_log10(labels=comma)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(color="orange", method="lm", se=TRUE)+
  xlab("Weights")+
  ylab("Density")

ggsave("Figures/Supplementary_figures/density_vs_weights_across_grids_15p.png", 
       width=6.8, height=6.5, units="in")

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

# Make a histogram of the densities of the 9,700 species
ggplot(weighted_dens, aes(x=10^density))+
  geom_histogram(bins=50, color="black", fill="lightblue")+
  scale_x_log10(labels=label_comma(accuracy=1))+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Density (birds per square mile)")+
  ylab("Number of species")

ggsave("Figures/Supplementary_figures/histogram_of_weighted_densities_15p.png",
       width=6.8, height=5.6, units="in")

# make a histogram of the standard errors, incorporarting imputation
ggplot(weighted_dens, aes(x=se_imp))+
  geom_histogram(bins=70, color="black", fill="lightblue")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Standard error of weighted density")+
  ylab("Number of species")

ggsave("Figures/Supplementary_figures/histogram_of_density_standard_errors_15p.png",
       width=6.8, height=5.6, units="in")

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
all_bird_total_results <- tibble(birds_non = spp_sum_non,
                                 birds_imp = spp_sum_imp) %>%
  mutate(mean_non=round(mean(spp_sum_non)),
         mean_imp=round(mean(spp_sum_imp)),
         median_non=round(median(spp_sum_non)),
         median_imp=round(median(spp_sum_imp))) %>%
  mutate(hdi_95_lower_non=round(hdi(spp_sum_non, 0.95))[[1]],
         hdi_95_upper_non=round(hdi(spp_sum_non, 0.95))[[2]],
         hdi_80_lower_non=round(hdi(spp_sum_non, 0.80))[[1]],
         hdi_80_upper_non=round(hdi(spp_sum_non, 0.80))[[2]],
         hdi_95_lower_imp=round(hdi(spp_sum_imp, 0.95))[[1]],
         hdi_95_upper_imp=round(hdi(spp_sum_imp, 0.95))[[2]],
         hdi_80_lower_imp=round(hdi(spp_sum_imp, 0.80))[[1]],
         hdi_80_upper_imp=round(hdi(spp_sum_imp, 0.80))[[2]]) %>%
  mutate(conf_95_lower_non=round((quantile(spp_sum_non, c(.025))))[[1]],
         conf_95_upper_non=round((quantile(spp_sum_non, c(.975))))[[1]],
         conf_80_lower_non=round((quantile(spp_sum_non, c(.10))))[[1]],
         conf_80_upper_non=round((quantile(spp_sum_non, c(.90))))[[1]],
         conf_95_lower_imp=round((quantile(spp_sum_imp, c(.025))))[[1]],
         conf_95_upper_imp=round((quantile(spp_sum_imp, c(.975))))[[1]],
         conf_80_lower_imp=round((quantile(spp_sum_imp, c(.10))))[[1]],
         conf_80_upper_imp=round((quantile(spp_sum_imp, c(.90))))[[1]])

saveRDS(all_bird_total_results, "Intermediate data files/imputation_summaries/reduced_15p/total_number_of_individual_birds_15p.RDS")


# want to write out the data to make the distributions for a 
# select number of species
# which can be presented in the main manuscript
# save each species
# as a separate RDS and then they can be joined later on
save_example_species_function <- function(example_species) {
  
  species_col_number <- sppdat %>%
    mutate(row_id=1:nrow(.)) %>%
    dplyr::filter(ebird_COMMON_NAME == example_species) %>%
    .$row_id
  
  species_dat <- tibble(non = spp_num_non[[species_col_number]],
                        imp = spp_num_imp[[species_col_number]]) %>%
    mutate(ebird_COMMON_NAME=sppdat$ebird_COMMON_NAME[[species_col_number]])
  
  saveRDS(species_dat, paste0("Intermediate data files/imputation_summaries/reduced_15p/", example_species, ".RDS"))
  
}

example_species <- c("Yellow-tailed Black-Cockatoo", "Ashy Prinia", "Midget Flowerpecker", 
                     "Yellow-vented Bulbul", "Coppersmith Barbet", "Osprey", "Acorn Woodpecker",
                     "Northern Wheatear", "Green Heron", "Ring-billed Gull", "Laughing Kookaburra",
                     "Eastern Wood-Pewee")

lapply(example_species, save_example_species_function)

# now need to do this but for each species
# and calculate the median as the abundance estimate
species_summary_function <- function(col_number){
  
  species_dat <- tibble(non = spp_num_non[[col_number]],
                        imp = spp_num_imp[[col_number]]) %>%
    mutate(ebird_COMMON_NAME=sppdat$ebird_COMMON_NAME[[col_number]])
  
  species_summary_df <- data.frame(ebird_COMMON_NAME=sppdat$ebird_COMMON_NAME[[col_number]],
                                   mean_non=round(mean(species_dat$non)),
                                   mean_imp=round(mean(species_dat$imp)),
                                   median_non=round(median(species_dat$non)),
                                   median_imp=round(median(species_dat$imp)),
                                   hdi_95_lower_non=round(hdi(species_dat$non, 0.95))[[1]],
                                   hdi_95_upper_non=round(hdi(species_dat$non, 0.95))[[2]],
                                   hdi_80_lower_non=round(hdi(species_dat$non, 0.80))[[1]],
                                   hdi_80_upper_non=round(hdi(species_dat$imp, 0.80))[[2]],
                                   hdi_95_lower_imp=round(hdi(species_dat$imp, 0.95))[[1]],
                                   hdi_95_upper_imp=round(hdi(species_dat$imp, 0.95))[[2]],
                                   hdi_80_lower_imp=round(hdi(species_dat$imp, 0.80))[[1]],
                                   hdi_80_upper_imp=round(hdi(species_dat$imp, 0.80))[[2]],
                                   conf_95_lower_non=round((quantile(species_dat$non, c(.025))))[[1]],
                                   conf_95_upper_non=round((quantile(species_dat$non, c(.975))))[[1]],
                                   conf_80_lower_non=round((quantile(species_dat$non, c(.10))))[[1]],
                                   conf_80_upper_non=round((quantile(species_dat$non, c(.90))))[[1]],
                                   conf_95_lower_imp=round((quantile(species_dat$imp, c(.025))))[[1]],
                                   conf_95_upper_imp=round((quantile(species_dat$imp, c(.975))))[[1]],
                                   conf_80_lower_imp=round((quantile(species_dat$imp, c(.10))))[[1]],
                                   conf_80_upper_imp=round((quantile(species_dat$imp, c(.90))))[[1]])
  
  return(species_summary_df)
}

# now get species summary as well
species_specific_results <- bind_rows(lapply(c(1:nrow(sppdat)), species_summary_function))

saveRDS(species_specific_results, "Intermediate data files/imputation_summaries/reduced_15p/species_specific_results_summary_15p.RDS")

# now get ORDER-level summaries for plotting
temp_dat <- sppdat %>%
  left_join(., clements) %>%
  group_by(ebird_COMMON_NAME) %>%
  slice(1)

order_summary_function <- function(order_name){
  
  temp_dat2 <- temp_dat %>%
    dplyr::filter(order==order_name)
  
  # abundance for all the species as ditribuitons of 10,000
  list_dens_non <- map2(temp_dat2$density, temp_dat2$se_non, ~ rlnorm(1e4, .x*log(10), .y*log(10))) %>% purrr::map(~.[order(.)]) 
  list_dens_imp <- map2(temp_dat2$density, temp_dat2$se_imp, ~ rlnorm(1e4, .x*log(10), .y*log(10))) %>% purrr::map(~.[order(.)])
  
  # 10,000 draws for all the spp 
  spp_num_non <- map2_dfc(list_dens_non, temp_dat2$total_area, ~.x*.y)              
  spp_num_imp <- map2_dfc(list_dens_imp, temp_dat2$total_area, ~.x*.y)  
  
  # summing these draws
  # summing is the step - means between two things become different
  spp_sum_non <- as.vector(rowSums(spp_num_non))
  spp_sum_imp <- as.vector(rowSums(spp_num_imp))
  
  
  order_results <- tibble(birds_non = spp_sum_non,
                          birds_imp = spp_sum_imp) %>%
    mutate(order=order_name) %>%
    mutate(mean_non=round(mean(spp_sum_non)),
           mean_imp=round(mean(spp_sum_imp)),
           median_non=round(median(spp_sum_non)),
           median_imp=round(median(spp_sum_imp))) %>%
    mutate(hdi_95_lower_non=round(hdi(spp_sum_non, 0.95))[[1]],
           hdi_95_upper_non=round(hdi(spp_sum_non, 0.95))[[2]],
           hdi_80_lower_non=round(hdi(spp_sum_non, 0.80))[[1]],
           hdi_80_upper_non=round(hdi(spp_sum_non, 0.80))[[2]],
           hdi_95_lower_imp=round(hdi(spp_sum_imp, 0.95))[[1]],
           hdi_95_upper_imp=round(hdi(spp_sum_imp, 0.95))[[2]],
           hdi_80_lower_imp=round(hdi(spp_sum_imp, 0.80))[[1]],
           hdi_80_upper_imp=round(hdi(spp_sum_imp, 0.80))[[2]]) %>%
    mutate(conf_95_lower_non=round((quantile(spp_sum_non, c(.025))))[[1]],
           conf_95_upper_non=round((quantile(spp_sum_non, c(.975))))[[1]],
           conf_80_lower_non=round((quantile(spp_sum_non, c(.10))))[[1]],
           conf_80_upper_non=round((quantile(spp_sum_non, c(.90))))[[1]],
           conf_95_lower_imp=round((quantile(spp_sum_imp, c(.025))))[[1]],
           conf_95_upper_imp=round((quantile(spp_sum_imp, c(.975))))[[1]],
           conf_80_lower_imp=round((quantile(spp_sum_imp, c(.10))))[[1]],
           conf_80_upper_imp=round((quantile(spp_sum_imp, c(.90))))[[1]])
  
}

# apply this function over each order in our data
order_results <- bind_rows(lapply(unique(temp_dat$order), order_summary_function))

saveRDS(order_results, "Intermediate data files/imputation_summaries/reduced_15p/order_specific_results_summary_15p.RDS")

# now get FAMILY-level summaries for plotting
temp_dat <- sppdat %>%
  left_join(., clements) %>%
  group_by(ebird_COMMON_NAME) %>%
  slice(1)

family_summary_function <- function(family_name){
  
  temp_dat2 <- temp_dat %>%
    dplyr::filter(family==family_name)
  
  # abundance for all the species as ditribuitons of 10,000
  list_dens_non <- map2(temp_dat2$density, temp_dat2$se_non, ~ rlnorm(1e4, .x*log(10), .y*log(10))) %>% purrr::map(~.[order(.)]) 
  list_dens_imp <- map2(temp_dat2$density, temp_dat2$se_imp, ~ rlnorm(1e4, .x*log(10), .y*log(10))) %>% purrr::map(~.[order(.)])
  
  # 10,000 draws for all the spp 
  spp_num_non <- map2_dfc(list_dens_non, temp_dat2$total_area, ~.x*.y)              
  spp_num_imp <- map2_dfc(list_dens_imp, temp_dat2$total_area, ~.x*.y)  
  
  # summing these draws
  # summing is the step - means between two things become different
  spp_sum_non <- as.vector(rowSums(spp_num_non))
  spp_sum_imp <- as.vector(rowSums(spp_num_imp))
  
  
  family_results <- tibble(birds_non = spp_sum_non,
                           birds_imp = spp_sum_imp) %>%
    mutate(family=family_name) %>%
    mutate(mean_non=round(mean(spp_sum_non)),
           mean_imp=round(mean(spp_sum_imp)),
           median_non=round(median(spp_sum_non)),
           median_imp=round(median(spp_sum_imp))) %>%
    mutate(hdi_95_lower_non=round(hdi(spp_sum_non, 0.95))[[1]],
           hdi_95_upper_non=round(hdi(spp_sum_non, 0.95))[[2]],
           hdi_80_lower_non=round(hdi(spp_sum_non, 0.80))[[1]],
           hdi_80_upper_non=round(hdi(spp_sum_non, 0.80))[[2]],
           hdi_95_lower_imp=round(hdi(spp_sum_imp, 0.95))[[1]],
           hdi_95_upper_imp=round(hdi(spp_sum_imp, 0.95))[[2]],
           hdi_80_lower_imp=round(hdi(spp_sum_imp, 0.80))[[1]],
           hdi_80_upper_imp=round(hdi(spp_sum_imp, 0.80))[[2]]) %>%
    mutate(conf_95_lower_non=round((quantile(spp_sum_non, c(.025))))[[1]],
           conf_95_upper_non=round((quantile(spp_sum_non, c(.975))))[[1]],
           conf_80_lower_non=round((quantile(spp_sum_non, c(.10))))[[1]],
           conf_80_upper_non=round((quantile(spp_sum_non, c(.90))))[[1]],
           conf_95_lower_imp=round((quantile(spp_sum_imp, c(.025))))[[1]],
           conf_95_upper_imp=round((quantile(spp_sum_imp, c(.975))))[[1]],
           conf_80_lower_imp=round((quantile(spp_sum_imp, c(.10))))[[1]],
           conf_80_upper_imp=round((quantile(spp_sum_imp, c(.90))))[[1]])
  
}

# apply this function over each order in our data
family_results <- bind_rows(lapply(unique(temp_dat$family), family_summary_function))

saveRDS(family_results, "Intermediate data files/imputation_summaries/reduced_15p/family_specific_results_summary_15p.RDS")



##### Get summaries at different 'organization' levels
##### similar to how order was done above

# now get biogeographic realm-level data
realm_dat <- sppdat %>%
  left_join(., clements) %>%
  group_by(ebird_COMMON_NAME) %>%
  slice(1) %>%
  left_join(., niches %>%
              dplyr::select(TipLabel, Realm)) %>%
  dplyr::filter(complete.cases(Realm))

realm_summary_function <- function(realm_name){
  
  dat <- realm_dat %>%
    dplyr::filter(Realm==realm_name)
  
  # abundance for all the species as ditribuitons of 10,000
  list_dens_non <- map2(dat$density, dat$se_non, ~ rlnorm(1e4, .x*log(10), .y*log(10))) %>% 
    purrr::map(~.[order(.)]) 
  list_dens_imp <- map2(dat$density, dat$se_imp, ~ rlnorm(1e4, .x*log(10), .y*log(10))) %>% 
    purrr::map(~.[order(.)])
  
  # 10,000 draws for all the spp 
  spp_num_non <- map2_dfc(list_dens_non, dat$total_area, ~.x*.y)              
  spp_num_imp <- map2_dfc(list_dens_imp, dat$total_area, ~.x*.y)  
  
  # summing these draws
  # summing is the step - means between two things become different
  spp_sum_non <- as.vector(rowSums(spp_num_non))
  spp_sum_imp <- as.vector(rowSums(spp_num_imp))
  
  
  realm_results <- tibble(birds_non = spp_sum_non,
                          birds_imp = spp_sum_imp) %>%
    mutate(realm=realm_name) %>%
    mutate(mean_non=round(mean(spp_sum_non)),
           mean_imp=round(mean(spp_sum_imp)),
           median_non=round(median(spp_sum_non)),
           median_imp=round(median(spp_sum_imp))) %>%
    mutate(hdi_95_lower_non=round(hdi(spp_sum_non, 0.95))[[1]],
           hdi_95_upper_non=round(hdi(spp_sum_non, 0.95))[[2]],
           hdi_80_lower_non=round(hdi(spp_sum_non, 0.80))[[1]],
           hdi_80_upper_non=round(hdi(spp_sum_non, 0.80))[[2]],
           hdi_95_lower_imp=round(hdi(spp_sum_imp, 0.95))[[1]],
           hdi_95_upper_imp=round(hdi(spp_sum_imp, 0.95))[[2]],
           hdi_80_lower_imp=round(hdi(spp_sum_imp, 0.80))[[1]],
           hdi_80_upper_imp=round(hdi(spp_sum_imp, 0.80))[[2]]) %>%
    mutate(conf_95_lower_non=round((quantile(spp_sum_non, c(.025))))[[1]],
           conf_95_upper_non=round((quantile(spp_sum_non, c(.975))))[[1]],
           conf_80_lower_non=round((quantile(spp_sum_non, c(.10))))[[1]],
           conf_80_upper_non=round((quantile(spp_sum_non, c(.90))))[[1]],
           conf_95_lower_imp=round((quantile(spp_sum_imp, c(.025))))[[1]],
           conf_95_upper_imp=round((quantile(spp_sum_imp, c(.975))))[[1]],
           conf_80_lower_imp=round((quantile(spp_sum_imp, c(.10))))[[1]],
           conf_80_upper_imp=round((quantile(spp_sum_imp, c(.90))))[[1]]) %>%
    mutate(number_of_species=nrow(dat))
  
}

# apply this function over each realm in our data
realm_results <- bind_rows(lapply(unique(realm_dat$Realm), realm_summary_function))

saveRDS(realm_results, "Intermediate data files/imputation_summaries/reduced_15p/realm_specific_results_summary_15p.RDS")


# repeat the above but for trophic level
trophic_level_dat <- sppdat %>%
  left_join(., clements) %>%
  group_by(ebird_COMMON_NAME) %>%
  slice(1) %>%
  left_join(., niches %>%
              dplyr::select(TipLabel, TrophicLevel)) %>%
  dplyr::filter(complete.cases(TrophicLevel))

trophic_level_summary_function <- function(trophic_level_name){
  
  dat <- trophic_level_dat %>%
    dplyr::filter(TrophicLevel==trophic_level_name)
  
  # abundance for all the species as ditribuitons of 10,000
  list_dens_non <- map2(dat$density, dat$se_non, ~ rlnorm(1e4, .x*log(10), .y*log(10))) %>% 
    purrr::map(~.[order(.)]) 
  list_dens_imp <- map2(dat$density, dat$se_imp, ~ rlnorm(1e4, .x*log(10), .y*log(10))) %>% 
    purrr::map(~.[order(.)])
  
  # 10,000 draws for all the spp 
  spp_num_non <- map2_dfc(list_dens_non, dat$total_area, ~.x*.y)              
  spp_num_imp <- map2_dfc(list_dens_imp, dat$total_area, ~.x*.y)  
  
  # summing these draws
  # summing is the step - means between two things become different
  spp_sum_non <- as.vector(rowSums(spp_num_non))
  spp_sum_imp <- as.vector(rowSums(spp_num_imp))
  
  
  realm_results <- tibble(birds_non = spp_sum_non,
                          birds_imp = spp_sum_imp) %>%
    mutate(trophic_level=trophic_level_name) %>%
    mutate(mean_non=round(mean(spp_sum_non)),
           mean_imp=round(mean(spp_sum_imp)),
           median_non=round(median(spp_sum_non)),
           median_imp=round(median(spp_sum_imp))) %>%
    mutate(hdi_95_lower_non=round(hdi(spp_sum_non, 0.95))[[1]],
           hdi_95_upper_non=round(hdi(spp_sum_non, 0.95))[[2]],
           hdi_80_lower_non=round(hdi(spp_sum_non, 0.80))[[1]],
           hdi_80_upper_non=round(hdi(spp_sum_non, 0.80))[[2]],
           hdi_95_lower_imp=round(hdi(spp_sum_imp, 0.95))[[1]],
           hdi_95_upper_imp=round(hdi(spp_sum_imp, 0.95))[[2]],
           hdi_80_lower_imp=round(hdi(spp_sum_imp, 0.80))[[1]],
           hdi_80_upper_imp=round(hdi(spp_sum_imp, 0.80))[[2]]) %>%
    mutate(conf_95_lower_non=round((quantile(spp_sum_non, c(.025))))[[1]],
           conf_95_upper_non=round((quantile(spp_sum_non, c(.975))))[[1]],
           conf_80_lower_non=round((quantile(spp_sum_non, c(.10))))[[1]],
           conf_80_upper_non=round((quantile(spp_sum_non, c(.90))))[[1]],
           conf_95_lower_imp=round((quantile(spp_sum_imp, c(.025))))[[1]],
           conf_95_upper_imp=round((quantile(spp_sum_imp, c(.975))))[[1]],
           conf_80_lower_imp=round((quantile(spp_sum_imp, c(.10))))[[1]],
           conf_80_upper_imp=round((quantile(spp_sum_imp, c(.90))))[[1]]) %>%
    mutate(number_of_species=nrow(dat))
  
}

# apply this function over each level in our data
trophic_level_results <- bind_rows(lapply(unique(trophic_level_dat$TrophicLevel), trophic_level_summary_function))

saveRDS(trophic_level_results, "Intermediate data files/imputation_summaries/reduced_15p/trophic_level_specific_results_summary_15p.RDS")

# repeat the above but for trophic niche
trophic_niche_dat <- sppdat %>%
  left_join(., clements) %>%
  group_by(ebird_COMMON_NAME) %>%
  slice(1) %>%
  left_join(., niches %>%
              dplyr::select(TipLabel, TrophicNiche)) %>%
  dplyr::filter(complete.cases(TrophicNiche))

trophic_niche_summary_function <- function(trophic_niche_name){
  
  dat <- trophic_niche_dat %>%
    dplyr::filter(TrophicNiche==trophic_niche_name)
  
  # abundance for all the species as ditribuitons of 10,000
  list_dens_non <- map2(dat$density, dat$se_non, ~ rlnorm(1e4, .x*log(10), .y*log(10))) %>% 
    purrr::map(~.[order(.)]) 
  list_dens_imp <- map2(dat$density, dat$se_imp, ~ rlnorm(1e4, .x*log(10), .y*log(10))) %>% 
    purrr::map(~.[order(.)])
  
  # 10,000 draws for all the spp 
  spp_num_non <- map2_dfc(list_dens_non, dat$total_area, ~.x*.y)              
  spp_num_imp <- map2_dfc(list_dens_imp, dat$total_area, ~.x*.y)  
  
  # summing these draws
  # summing is the step - means between two things become different
  spp_sum_non <- as.vector(rowSums(spp_num_non))
  spp_sum_imp <- as.vector(rowSums(spp_num_imp))
  
  
  realm_results <- tibble(birds_non = spp_sum_non,
                          birds_imp = spp_sum_imp) %>%
    mutate(trophic_niche=trophic_niche_name) %>%
    mutate(mean_non=round(mean(spp_sum_non)),
           mean_imp=round(mean(spp_sum_imp)),
           median_non=round(median(spp_sum_non)),
           median_imp=round(median(spp_sum_imp))) %>%
    mutate(hdi_95_lower_non=round(hdi(spp_sum_non, 0.95))[[1]],
           hdi_95_upper_non=round(hdi(spp_sum_non, 0.95))[[2]],
           hdi_80_lower_non=round(hdi(spp_sum_non, 0.80))[[1]],
           hdi_80_upper_non=round(hdi(spp_sum_non, 0.80))[[2]],
           hdi_95_lower_imp=round(hdi(spp_sum_imp, 0.95))[[1]],
           hdi_95_upper_imp=round(hdi(spp_sum_imp, 0.95))[[2]],
           hdi_80_lower_imp=round(hdi(spp_sum_imp, 0.80))[[1]],
           hdi_80_upper_imp=round(hdi(spp_sum_imp, 0.80))[[2]]) %>%
    mutate(conf_95_lower_non=round((quantile(spp_sum_non, c(.025))))[[1]],
           conf_95_upper_non=round((quantile(spp_sum_non, c(.975))))[[1]],
           conf_80_lower_non=round((quantile(spp_sum_non, c(.10))))[[1]],
           conf_80_upper_non=round((quantile(spp_sum_non, c(.90))))[[1]],
           conf_95_lower_imp=round((quantile(spp_sum_imp, c(.025))))[[1]],
           conf_95_upper_imp=round((quantile(spp_sum_imp, c(.975))))[[1]],
           conf_80_lower_imp=round((quantile(spp_sum_imp, c(.10))))[[1]],
           conf_80_upper_imp=round((quantile(spp_sum_imp, c(.90))))[[1]]) %>%
    mutate(number_of_species=nrow(dat))
  
}

# apply this function over each level in our data
trophic_niche_results <- bind_rows(lapply(unique(trophic_niche_dat$TrophicNiche), trophic_niche_summary_function))

saveRDS(trophic_niche_results, "Intermediate data files/imputation_summaries/reduced_15p/trophic_niche_specific_results_summary_15p.RDS")

# repeat the above but for foraging niche
foraging_niche_dat <- sppdat %>%
  left_join(., clements) %>%
  group_by(ebird_COMMON_NAME) %>%
  slice(1) %>%
  left_join(., niches %>%
              dplyr::select(TipLabel, ForagingNiche)) %>%
  dplyr::filter(complete.cases(ForagingNiche))

foraging_niche_summary_function <- function(foraging_niche_name){
  
  dat <- foraging_niche_dat %>%
    dplyr::filter(ForagingNiche==foraging_niche_name)
  
  # abundance for all the species as ditribuitons of 10,000
  list_dens_non <- map2(dat$density, dat$se_non, ~ rlnorm(1e4, .x*log(10), .y*log(10))) %>% 
    purrr::map(~.[order(.)]) 
  list_dens_imp <- map2(dat$density, dat$se_imp, ~ rlnorm(1e4, .x*log(10), .y*log(10))) %>% 
    purrr::map(~.[order(.)])
  
  # 10,000 draws for all the spp 
  spp_num_non <- map2_dfc(list_dens_non, dat$total_area, ~.x*.y)              
  spp_num_imp <- map2_dfc(list_dens_imp, dat$total_area, ~.x*.y)  
  
  # summing these draws
  # summing is the step - means between two things become different
  spp_sum_non <- as.vector(rowSums(spp_num_non))
  spp_sum_imp <- as.vector(rowSums(spp_num_imp))
  
  
  realm_results <- tibble(birds_non = spp_sum_non,
                          birds_imp = spp_sum_imp) %>%
    mutate(foraging_niche=foraging_niche_name) %>%
    mutate(mean_non=round(mean(spp_sum_non)),
           mean_imp=round(mean(spp_sum_imp)),
           median_non=round(median(spp_sum_non)),
           median_imp=round(median(spp_sum_imp))) %>%
    mutate(hdi_95_lower_non=round(hdi(spp_sum_non, 0.95))[[1]],
           hdi_95_upper_non=round(hdi(spp_sum_non, 0.95))[[2]],
           hdi_80_lower_non=round(hdi(spp_sum_non, 0.80))[[1]],
           hdi_80_upper_non=round(hdi(spp_sum_non, 0.80))[[2]],
           hdi_95_lower_imp=round(hdi(spp_sum_imp, 0.95))[[1]],
           hdi_95_upper_imp=round(hdi(spp_sum_imp, 0.95))[[2]],
           hdi_80_lower_imp=round(hdi(spp_sum_imp, 0.80))[[1]],
           hdi_80_upper_imp=round(hdi(spp_sum_imp, 0.80))[[2]]) %>%
    mutate(conf_95_lower_non=round((quantile(spp_sum_non, c(.025))))[[1]],
           conf_95_upper_non=round((quantile(spp_sum_non, c(.975))))[[1]],
           conf_80_lower_non=round((quantile(spp_sum_non, c(.10))))[[1]],
           conf_80_upper_non=round((quantile(spp_sum_non, c(.90))))[[1]],
           conf_95_lower_imp=round((quantile(spp_sum_imp, c(.025))))[[1]],
           conf_95_upper_imp=round((quantile(spp_sum_imp, c(.975))))[[1]],
           conf_80_lower_imp=round((quantile(spp_sum_imp, c(.10))))[[1]],
           conf_80_upper_imp=round((quantile(spp_sum_imp, c(.90))))[[1]]) %>%
    mutate(number_of_species=nrow(dat))
  
}

# apply this function over each level in our data
foraging_niche_results <- bind_rows(lapply(unique(foraging_niche_dat$ForagingNiche), foraging_niche_summary_function))

saveRDS(foraging_niche_results, "Intermediate data files/imputation_summaries/reduced_15p/foraging_niche_specific_summary_results_15p.RDS")


# MAke a "Table" That will hold all the species-specific results for publication
table_s1 <- species_specific_results %>%
  dplyr::select(1, 5, 18, 19) %>%
  left_join(., sppdat %>%
              dplyr::select(ebird_COMMON_NAME, range_adjusted)) %>%
  left_join(., dat %>%
              dplyr::filter(complete.cases(density)) %>%
              dplyr::select(ebird_COMMON_NAME) %>%
              distinct() %>%
              mutate(Training_species="TRUE")) %>%
  replace_na(list(Training_species="FALSE")) %>%
  left_join(., clements %>%
              dplyr::select(1:4) %>%
              distinct()) %>%
  dplyr::select(1, 7:9, 2:6) %>%
  arrange(desc(median_imp))

write_csv(table_s1, "Tables/all_species_summary_table_15p.csv")

final_data <- read.csv("Tables/all_species_summary_table_15p.csv")
