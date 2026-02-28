## This is a script to prepare the various data
## from training datasets for modelling
## it also looks at the relationships
## between modelled abundances predicted from a GLM
## and more simple relative abundance measures taken by summing total birds and
## total effort, etc.
## there are three types of 'training data' and each of those are treated in turn here
## and then combined into a single dataset

# packages
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggcorrplot)
library(corrplot)
library(scales)

# read in states data
states <- read_csv("Data/states.csv")

# split for gb and combine among the multiple states
# because we don't have this data available to us
# we only have all of Great Britain estimates
gb_states <- states %>%
  dplyr::filter(STATE_CODE %in% c("GB-ENG", "GB-SCT", "GB-WLS", "GB-NIR")) %>%
  mutate(STATE_CODE2="GB") %>%
  mutate(area_square_miles2=sum(area_square_miles)) %>%
  dplyr::select(STATE_CODE2, area_square_miles2) %>%
  rename(STATE_CODE=STATE_CODE2) %>%
  rename(area_square_miles=area_square_miles2) %>%
  distinct()

# set working directory to the analyzed eBird abundances
setwd("Data/analyzed_eBird_abundances")

# read in data
# first GB separately
# and start to clean the data getting rid of some
# of the modelled estimates that are too wacky
# related to model size
GB_df <- data.frame(files=list.files()) %>%
  mutate(files=as.character(files)) %>%
  mutate(GB=startsWith(.$files, "GB")) %>%
  dplyr::filter(GB=="TRUE") %>%
  .$files %>%
  map_dfr(readRDS) %>%
  bind_rows() %>%
  mutate(abund_per_time=sum_of_abundance/total_effort_time) %>%
  # only look at models that converged
  dplyr::filter(model_converged=="TRUE") %>%
  # still lots of noise in the modelled abundances
  # and now filter out any abundance estimates that are > .99 quantile
  # and < .01 quantile. Simply 'cleaning' the noisy possible data points
  mutate(quantile.99=quantile(model_fit, 0.99)) %>%
  mutate(quantile.01=quantile(model_fit, 0.01)) %>%
  mutate(within_filtered_realm=ifelse(model_fit < quantile.99 & model_fit > quantile.01, "TRUE", "FALSE")) %>%
  dplyr::filter(within_filtered_realm=="TRUE")

# Now do the same thing for North American data
NA_df <- data.frame(files=list.files()) %>%
  mutate(files=as.character(files)) %>%
  mutate(North_America=startsWith(.$files, "NA")) %>%
  dplyr::filter(North_America=="TRUE") %>%
  .$files %>%
  map_dfr(readRDS) %>%
  bind_rows() %>%
  mutate(abund_per_time=sum_of_abundance/total_effort_time) %>%
  # only look at models that converged
  dplyr::filter(model_converged=="TRUE") %>%
  # still lots of noise in the modelled abundances
  # and now filter out any abundance estimates that are > .99 quantile
  # and < .01 quantile. Simply 'cleaning' the noisy possible data points
  mutate(quantile.99=quantile(model_fit, 0.99)) %>%
  mutate(quantile.01=quantile(model_fit, 0.01)) %>%
  mutate(within_filtered_realm=ifelse(model_fit < quantile.99 & model_fit > quantile.01, "TRUE", "FALSE")) %>%
  dplyr::filter(within_filtered_realm=="TRUE")

# Now do the same thing for the Birdlife data
Birdlife_df <- data.frame(files=list.files()) %>%
  mutate(files=as.character(files)) %>%
  mutate(Birdlife=startsWith(.$files, "Birdlife")) %>%
  dplyr::filter(Birdlife=="TRUE") %>%
  .$files %>%
  map_dfr(readRDS) %>%
  bind_rows() %>%
  mutate(abund_per_time=sum_of_abundance/total_effort_time) %>%
  # only look at models that converged
  dplyr::filter(model_converged=="TRUE") %>%
  # still lots of noise in the modelled abundances
  # and now filter out any abundance estimates that are > .99 quantile
  # and < .01 quantile. Simply 'cleaning' the noisy possible data points
  mutate(quantile.99=quantile(model_fit, 0.99)) %>%
  mutate(quantile.01=quantile(model_fit, 0.01)) %>%
  mutate(within_filtered_realm=ifelse(model_fit < quantile.99 & model_fit > quantile.01, "TRUE", "FALSE")) %>%
  dplyr::filter(within_filtered_realm=="TRUE") %>%
  mutate(Estimate=as.numeric(as.character(Estimate))) %>%
  mutate(Lower=as.numeric(as.character(Lower))) %>%
  mutate(Upper=as.numeric(as.character(Upper)))

setwd('..')
setwd('..')

# preliminary investigation into the correlation among our 
# potential different measures of relative abundance
# as calculated from eBird
## make a plot to look at the correlation between our three different 
## potential measures of relative abundance
## in theory it *shouldn't* matter as we are simply training the eBird sampling approach
## on some external data
mat <- bind_rows(GB_df, NA_df, Birdlife_df) %>%
  ungroup() %>%
  dplyr::select(mean_of_abundance, model_fit, abund_per_time) %>%
  rename(modelled_abundance=model_fit) %>%
  rename(abundance_per_time=abund_per_time) %>%
  as.matrix() %>%
  cor(.)

ggcorrplot(mat,
           ggtheme = ggplot2::theme_bw,
           lab=TRUE)

# now filter and clean the data for North America
# we filter only to summer months because all estimates are based on 
# breeding bird survey estimates

NA_mod_populations <- NA_df %>%
  dplyr::filter(MONTH %in% c("May", "Jun", "Jul", "Aug")) %>% 
  dplyr::select(COMMON_NAME, Population_Estimate, `Lower_95%_bound`, `Upper_95%_bound`, 
                STATE_CODE, STATE_CODE_BCR) %>% 
  distinct() %>% 
  group_by(COMMON_NAME, STATE_CODE) %>% 
  summarize(Population_Estimate=sum(Population_Estimate),
            `Lower_95%_bound`=sum(`Lower_95%_bound`),
            `Upper_95%_bound`=sum(`Upper_95%_bound`)) %>%
  mutate(pop_var = (`Upper_95%_bound` - `Lower_95%_bound`) / (2*qnorm(0.975))^2,
         pop_sd = sqrt(pop_var)) %>%
  rename(pop_up = `Upper_95%_bound`,
         pop_down = `Lower_95%_bound`) 


NA_mod_df <- NA_df %>%
  dplyr::filter(MONTH %in% c("May" , "Jun", "Jul", "Aug")) %>%
  # add a simple measure of abundance per time
  mutate(abund_per_time = sum_of_abundance/total_effort_time) %>%
  group_by(COMMON_NAME, STATE_CODE) %>%
  summarise(count_sd = n(),
            total_obs = sum(number_obs_in_model),
            mean_abundance_from_model = mean(model_fit, na.rm = TRUE),
            mean_of_abund = mean(mean_of_abundance),
            mean_abund_per_time = mean(abund_per_time),
            ln_sd_of_abund = log(sd(mean_of_abundance)) + 1/(2*(count_sd - 1)), # sd (time x space) of estaimtes of mean estiamtes with small sample size corrections (small sample size tend to underestiamte SD)
            sd_of_abund = exp(ln_sd_of_abund)) %>%
  left_join(., NA_mod_populations) %>%
  # join to state areas and calculate density
  left_join(., states, by="STATE_CODE") %>%
  mutate(density = Population_Estimate / area_square_miles,
    dens_up = pop_up / area_square_miles,
    dens_down = pop_down / area_square_miles,
    dens_sd = (dens_up - dens_down) / (2*qnorm(0.975)), # SD for density (this is actually SE = SD of estimates)
    #dens_sd2 = (log(dens_up) - log(dens_down)) / (2*qnorm(0.975)*exp((log(dens_up) + log(dens_down))/2)), # this does not work as some desnsity estimates are 0 - then get Inf
    dens_mid = (dens_up + dens_down)/2,
    species_state = paste(COMMON_NAME, STATE_CODE)) %>%
  mutate(SEASON="Breeding") %>% 
  dplyr::select(-state)

# More-or-less repeat the same thing, but for GB
# will split the 'breeding' and 'wintering' estimates
# and then combine them back together
# so some species can have both wintering and breeding estimates
GB_mod_df_breeding <- GB_df %>%
  mutate(STATE_CODE = "GB") %>%
  mutate(pop_var = (Upper - Lower) / (2*qnorm(0.975))^2,
         pop_mid = (Upper + Lower)/2,
         log_pop = log10(Estimate),
         log_pop_mid = (log10(Upper)+log10(Lower))/2) %>%
  dplyr::filter(MONTH %in% c("May" , "Jun", "Jul", "Aug")) %>%
  dplyr::filter(Season == "B") %>%
  # add a simple measure of abundance per time
  mutate(abund_per_time = sum_of_abundance/total_effort_time) %>%
  group_by(COMMON_NAME, STATE_CODE) %>%
  summarise(count_sd = n(),
            total_obs = sum(number_obs_in_model),
            mean_abundance_from_model = mean(model_fit, na.rm = TRUE),
            mean_of_abund = mean(mean_of_abundance),
            mean_abund_per_time = mean(abund_per_time),
            ln_sd_of_abund = log(sd(mean_of_abundance)) + 1/(2*(count_sd - 1)), # sd (time x space) of estaimtes of mean estiamtes with small sample size corrections (small sample size tend to underestiamte SD)
            sd_of_abund = exp(ln_sd_of_abund), # original scale
            Population_Estimate = mean(Estimate),
            pop_var = mean(pop_var), # we should take mean at var not SD
            pop_sd = sqrt(pop_var),
            pop_up = mean(Upper),
            pop_down = mean(Lower)) %>%
  # join to state areas and calculate density
  left_join(., gb_states, by="STATE_CODE") %>%
  mutate(density = Population_Estimate / area_square_miles,
         dens_up = pop_up / area_square_miles,
         dens_down = pop_down / area_square_miles,
         dens_sd = (dens_up - dens_down) / (2*qnorm(0.975)), # SD for density (this is actually SE = SD of estimates)
         #dens_sd2 = (log(dens_up) - log(dens_down)) / (2*qnorm(0.975)*exp((log(dens_up) + log(dens_down))/2)), # this does not work as some desnsity estimates are 0 - then get Inf
         dens_mid = (dens_up + dens_down)/2,
         species_state = paste(COMMON_NAME, STATE_CODE)) %>%
  mutate(SEASON="Breeding") %>%
  dplyr::select(colnames(NA_mod_df))

# for the wintering estimates
GB_mod_df_wintering <- GB_df %>%
  mutate(STATE_CODE = "GB") %>%
  mutate(pop_var = (Upper - Lower) / (2*qnorm(0.975))^2,
         pop_mid = (Upper + Lower)/2,
         log_pop = log10(Estimate),
         log_pop_mid = (log10(Upper)+log10(Lower))/2) %>%
  dplyr::filter(MONTH %in% c("Oct" , "Nov", "Dec", "Jan", "Feb")) %>%
  dplyr::filter(Season == "W") %>%
  # add a simple measure of abundance per time
  mutate(abund_per_time = sum_of_abundance/total_effort_time) %>%
  group_by(COMMON_NAME, STATE_CODE) %>%
  summarise(count_sd = n(),
            total_obs = sum(number_obs_in_model),
            mean_abundance_from_model = mean(model_fit, na.rm = TRUE),
            mean_of_abund = mean(mean_of_abundance),
            mean_abund_per_time = mean(abund_per_time),
            ln_sd_of_abund = log(sd(mean_of_abundance)) + 1/(2*(count_sd - 1)), # sd (time x space) of estaimtes of mean estiamtes with small sample size corrections (small sample size tend to underestiamte SD)
            sd_of_abund = exp(ln_sd_of_abund), # original scale
            Population_Estimate = mean(Estimate),
            pop_var = mean(pop_var), # we should take mean at var not SD
            pop_sd = sqrt(pop_var),
            pop_up = mean(Upper),
            pop_down = mean(Lower)) %>%
  # join to state areas and calculate density
  left_join(., gb_states, by="STATE_CODE") %>%
  mutate(density = Population_Estimate / area_square_miles,
         dens_up = pop_up / area_square_miles,
         dens_down = pop_down / area_square_miles,
         dens_sd = (dens_up - dens_down) / (2*qnorm(0.975)), # SD for density (this is actually SE = SD of estimates)
         #dens_sd2 = (log(dens_up) - log(dens_down)) / (2*qnorm(0.975)*exp((log(dens_up) + log(dens_down))/2)), # this does not work as some desnsity estimates are 0 - then get Inf
         dens_mid = (dens_up + dens_down)/2,
         species_state = paste(COMMON_NAME, STATE_CODE)) %>%
  mutate(SEASON="Wintering") %>%
  dplyr::select(colnames(NA_mod_df))

GB_mod_df <- bind_rows(GB_mod_df_breeding, GB_mod_df_wintering)

## Now pull together a dataframe for
## the Birdlife data
Birdlife_mod_df <- Birdlife_df %>%
  mutate(STATE_CODE = "Birdlife") %>%
  mutate(pop_var = (Upper - Lower) / (2*qnorm(0.975))^2,
         pop_mid = (Upper + Lower)/2,
         log_pop = log10(Estimate),
         log_pop_mid = (log10(Upper)+log10(Lower))/2) %>%
  # add a simple measure of abundance per time
  mutate(abund_per_time = sum_of_abundance/total_effort_time) %>%
  group_by(COMMON_NAME, STATE_CODE) %>%
  summarise(count_sd = n(),
            total_obs = sum(number_obs_in_model),
            mean_abundance_from_model = mean(model_fit, na.rm = TRUE),
            mean_of_abund = mean(mean_of_abundance),
            mean_abund_per_time = mean(abund_per_time),
            ln_sd_of_abund = log(sd(mean_of_abundance)) + 1/(2*(count_sd - 1)), # sd (time x space) of estaimtes of mean estiamtes with small sample size corrections (small sample size tend to underestiamte SD)
            sd_of_abund = exp(ln_sd_of_abund), # original scale
            Population_Estimate = mean(Estimate),
            pop_var = mean(pop_var), # we should take mean at var not SD
            pop_sd = sqrt(pop_var),
            pop_up = mean(Upper),
            pop_down = mean(Lower)) %>%
  # join to state areas and calculate density
  left_join(., distinct(dplyr::select(Birdlife_df, COMMON_NAME, area_square_miles)), by="COMMON_NAME") %>%
  mutate(density = Population_Estimate / area_square_miles,
         dens_up = pop_up / area_square_miles,
         dens_down = pop_down / area_square_miles,
         dens_sd = (dens_up - dens_down) / (2*qnorm(0.975)), # SD for density (this is actually SE = SD of estimates)
         #dens_sd2 = (log(dens_up) - log(dens_down)) / (2*qnorm(0.975)*exp((log(dens_up) + log(dens_down))/2)), # this does not work as some desnsity estimates are 0 - then get Inf
         dens_mid = (dens_up + dens_down)/2,
         species_state = paste(COMMON_NAME, STATE_CODE)) %>%
  mutate(SEASON="Year-round") %>%
  dplyr::select(colnames(NA_mod_df))



# this is another look at the correlation among our variables
# but this time it looks at the seasons with which the data were generated
# and takes the mean of these values as well
# so makes a bit more sense
# and is more correlated
mat2 <- bind_rows(NA_mod_df, GB_mod_df, Birdlife_mod_df) %>%
  ungroup() %>%
  dplyr::select(mean_of_abund, mean_abundance_from_model, mean_abund_per_time) %>%
  rename(`Abundance per time`=mean_abund_per_time) %>%
  rename(`Modelled abundance`=mean_abundance_from_model) %>%
  rename(`Mean abundance`=mean_of_abund) %>%
  as.matrix() %>%
  cor(.)

ggcorrplot(mat2,
           ggtheme = ggplot2::theme_bw,
           lab=TRUE)+
  theme(axis.text=element_text(color="black"))

ggsave("Figures/Supplementary_figures/correlation_matrix_of_abundance_measures_mean.png", 
       width=7, height=6.5, units="in")

# we can adjust 0 SD values - SD = 0 is coming from mean = 0 values (probably these will be elimated)
#match(which(mod_df$mean_of_abund == 0), which(mod_df$sd_of_abund == 0))
#convert SD on original to SD on log scale using delta method - use log_10 not ln

mod_df <- bind_rows(NA_mod_df, GB_mod_df, Birdlife_mod_df) %>%
  mutate(species = COMMON_NAME, #just I wanted to make this
  density_log10 = log10(density),
  density_SD_log10 = dens_sd/(log(10)*density), # delta method (this one is already SE = SD)
  abundance_log10 = log10(mean_of_abund),
  abundance_SD_log10 = (sd_of_abund/sqrt(count_sd))/(log(10)*mean_of_abund)) # delta method

#### CUTOFF of observations ####

obs_percentage <- function(percentage){
  spec_obs <- bind_rows(NA_mod_df, GB_mod_df, Birdlife_mod_df) %>% 
    rename(species = COMMON_NAME) %>%
    group_by(species) %>% 
    summarise(total_obs_sum = sum(total_obs, na.rm = TRUE)) %>%
    ungroup() %>% 
    filter(total_obs_sum >= quantile(total_obs_sum, percentage))
  
  mod_df <- mod_df %>% 
    filter(species %in% spec_obs$species)
  
  return(mod_df)
}

modf_2p <- obs_percentage(0.02)
modf_5p <- obs_percentage(0.05)
modf_10p <- obs_percentage(0.1)
modf_15p <- obs_percentage(0.15)

saveRDS(mod_df, "Intermediate data files/data_for_detectability_modelling.RDS")
saveRDS(modf_2p, "Intermediate data files/data_for_detectability_modelling_2p.RDS")
saveRDS(modf_5p, "Intermediate data files/data_for_detectability_modelling_5p.RDS")
saveRDS(modf_10p, "Intermediate data files/data_for_detectability_modelling_10p.RDS")
saveRDS(modf_15p, "Intermediate data files/data_for_detectability_modelling_15p.RDS")


#new version of function where I also get list of removed species

obs_percentage <- function(percentage){
  
  spec_obs <- bind_rows(NA_mod_df, GB_mod_df, Birdlife_mod_df) %>% 
    rename(species = COMMON_NAME) %>%
    group_by(species) %>% 
    summarise(total_obs_sum = sum(total_obs, na.rm = TRUE)) %>%
    ungroup()
  
  # cutoff value
  cutoff <- quantile(spec_obs$total_obs_sum, percentage)
  
  # species to keep
  keep_species <- spec_obs %>% 
    filter(total_obs_sum >= cutoff)
  
  # species removed
  removed_species <- spec_obs %>% 
    filter(total_obs_sum < cutoff)
  
  # filtered dataset
  filtered_df <- mod_df %>% 
    filter(species %in% keep_species$species)
  
  return(list(
    data = filtered_df,
    removed = removed_species
  ))
}

result <- obs_percentage(0.15)

filtered_data <- result$data
removed_species <- result$removed





# make some summary plots of these training data
# some of which may be useful for supplementary figures

# plot of the histogram of the areas used for calculation
mod_df %>%
  dplyr::select(species_state, area_square_miles) %>%
  distinct() %>%
  ggplot(., aes(x=area_square_miles))+
  geom_histogram(bins=50, color="black", fill="lightblue")+
  scale_x_log10(labels=comma)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Area (square miles)")+
  ylab("Number of observations")

ggsave("Figures/Supplementary_figures/area_of_range_used_to_estimate_populations.png",
       width=6.8, height=5.6, units="in")

# plot of overall mean abundance
# versus total abundance
# divided by total total effort
ggplot(mod_df, aes(mean_of_abund, mean_abund_per_time))+
  geom_point(color="black")+
  scale_x_log10(labels=comma)+
  scale_y_log10(labels=comma)+
  theme_bw()+
  geom_smooth(method="lm", color="orange")+
  theme(axis.text=element_text(color="black"))+
  xlab("Overall mean abundance per time (log10)")+
  ylab("Abundance per effort (log10)")

ggsave("Figures/Supplementary_figures/mean_abund_vs_abund_per_time.png",
       width=6.8, height=5.6, units="in")

# now plot the underlying data for a handful of species
# as an example plot, illustrating what is being done
# for the training data
mod_df %>%
  dplyr::filter(
    COMMON_NAME %in% c(
      "Blue-gray Gnatcatcher",
      "American Crow",
      "Black-capped Chickadee",
      "American Robin"
    )
  ) %>%
  ggplot(., aes(x = mean_of_abund, y = density, col = COMMON_NAME)) +
  geom_point() +
  geom_errorbar(mapping = aes(x = mean_of_abund, ymin = dens_down, ymax =
                                dens_up)) +
  scale_y_log10(labels=comma) +
  scale_x_log10(labels=comma) +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  facet_wrap( ~ COMMON_NAME, scales = "free")+
  xlab("eBird relative abundance (log10)")+
  ylab("Observed density (log10)")+
  ggtitle("normal abundance") +
  guides(color=FALSE)+
  theme(axis.text=element_text(color="black"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_brewer(palette = "Set2")


modf_2p %>%
  dplyr::filter(
    COMMON_NAME %in% c(
      "Blue-gray Gnatcatcher",
      "American Crow",
      "Black-capped Chickadee",
      "American Robin"
    )
  ) %>%
  ggplot(., aes(x = mean_of_abund, y = density, col = COMMON_NAME)) +
  geom_point() +
  geom_errorbar(mapping = aes(x = mean_of_abund, ymin = dens_down, ymax =
                                dens_up)) +
  scale_y_log10(labels=comma) +
  scale_x_log10(labels=comma) +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_bw() +
  facet_wrap( ~ COMMON_NAME, scales = "free")+
  xlab("eBird relative abundance (log10)")+
  ylab("Observed density (log10)")+
  ggtitle("Abundance cutoff at 15%") +
  guides(color=FALSE)+
  theme(axis.text=element_text(color="black"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_brewer(palette = "Set2")


ggsave("Figures/Supplementary_figures/mean_abund_vs_observed_density_example_2p.png",
       width=6.8, height=6.5, units="in")




####### TRYING STUFF ########

max(mod_df$mean_of_abund)



table_before <- table(mod_df$mean_of_abund)
table_after  <- table(mod_df_2p$mean_of_abund)

# Species that became singletons
became_singleton <- names(table_after[table_after == 1 & table_before[names(table_after)] > 1])

length(became_singleton)

library(lattice)
xtabs( ~ mod_df$STATE_CODE)

min(mod_df$total_obs)
nrow(mod_df[mod_df$total_obs == "1" ,])
