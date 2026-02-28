####################################
# prepare_file_for_imputation
####################################
## This is an intermediate R script
## that pulls the slopes calculated from Shinichi's "density_abundance_model.Rmd"
## file and then the color and body size data necessary for imputation

# packages
pacman::p_load(tidyverse,
               magrittr,
               dplyr,
               readr,
               broom,
               lme4,
               GGally,
               corrplot,               
               here,
               rapportools,
               mvtnorm,
               sf,
               ggcorrplot
)


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

# read in results of brms model
detectability <- read_csv("Intermediate data files/spp_all.csv") %>%
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

# remove large files above
rm(temporal_summary)
rm(combined_dat)



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

# read in grid cell land mass estimates
grids_land_area <- readRDS("Data/Spatial_data/grid_percent_land.RDS") %>%
  dplyr::select(ID, percentage) %>%
  dplyr::rename(grid_id=ID) %>%
  st_set_geometry(NULL)

joined_dat_final <- bind_rows(training_dat.2, data_for_imputation.2) %>%
  left_join(., grids_land_area)

saveRDS(joined_dat_final, "Intermediate data files/file_for_imputation2.RDS")

# to save time and not run all of the above
#joined_dat_final <- readRDS("Intermediate data files/file_for_imputation2.RDS")

# summarize the information
# for these relationships
# with intercept and slope
# and the missingness (at the species level) for each trait
# prepare data for imputation
# by logging or standardizing values
temp_dat <- joined_dat_final %>%
  ungroup() %>%
  dplyr::select(10:18) %>%
  mutate(max_distance_log10=log10(max_distance)) %>%
  mutate(max_brightness_log10=log10(max_brightness)) %>%
  dplyr::filter(complete.cases(intercept)) %>%
  distinct() %>%
  dplyr::select(-max_distance, -max_brightness)

ggpairs(temp_dat)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  theme(strip.text=element_text(size=5))

ggsave("Figures/Supplementary_figures/intercept_and_slope_trait_relationships.png", 
       height=7.8, width=9.1, units="in")

species_level <- joined_dat_final %>%
  ungroup() %>%
  dplyr::select(2, 10:14) %>%
  distinct()

# percent missingness for each column
sum(is.na(species_level$max_distance))/nrow(species_level)*100

sum(is.na(species_level$max_brightness))/nrow(species_level)*100

sum(is.na(species_level$mass_log10))/nrow(species_level)*100

sum(is.na(species_level$flock_size_log10))/nrow(species_level)*100

sum(is.na(species_level$IUCN_ordinal))/nrow(species_level)*100

# percent missingness of density at the species level
(1-joined_dat_final %>%
    ungroup() %>%
    dplyr::select(ebird_COMMON_NAME, density) %>%
    dplyr::filter(complete.cases(density)) %>%
    dplyr::select(ebird_COMMON_NAME) %>%
    distinct() %>%
    nrow()/nrow(species_level))*100

# percent missingness of density throughout the whole imputation
(1-joined_dat_final %>%
    ungroup() %>%
    dplyr::select(ebird_COMMON_NAME, density) %>%
    dplyr::filter(complete.cases(density)) %>%
    nrow()/nrow(joined_dat_final))*100

# now look at pairwise relationships between life history traits and
# rel abund and density
joined_dat_final %>%
  ungroup() %>%
  dplyr::select(10:14, 20:22) %>%
  dplyr::rename(density_se=se) %>%
  mutate(max_distance_log10=log10(max_distance)) %>%
  mutate(max_brightness_log10=log10(max_brightness)) %>%
  dplyr::select(-max_distance, -max_brightness) %>%
  dplyr::filter(complete.cases(log10_density)) %>%
  as.matrix() %>%
  cor(., use="pairwise.complete.obs") %>%
  ggcorrplot(., ggtheme=ggplot2::theme_bw, lab=TRUE)+
  theme(axis.text=element_text(color="black"))

ggsave("Figures/Supplementary_figures/density_and_rel_abund_trait_relationships.png", 
       height=6.5, width=6.5, units="in")
