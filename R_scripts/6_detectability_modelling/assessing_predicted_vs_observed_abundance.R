### This is an R script to test our detectability function
### with the training data
### and re-predict the estimated abundances/densities to assess
### against the originally observed estimates


# packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(readr)

# read in results of MCMCglmm model
detectability <- read_csv("Intermediate data files/spp_all.csv") %>%
  pivot_wider(names_from=type, values_from=c(estimate, se)) %>%
  rename(ebird_COMMON_NAME=species) %>%
  rename(intercept=estimate_intercept) %>%
  rename(slope=estimate_slope) %>%
  rename(intercept_se=se_intercept) %>%
  rename(slope_se=se_slope)


# read in original observed data used for detectability modelling
mod_df <- readRDS("Intermediate data files/data_for_detectability_modelling.RDS")

# checking for infinity
which(mod_df$density_log10 == -Inf)
m_inf <- which(mod_df$abundance_log10 == -Inf)
min(mod_df$abundance_log10[-m_inf]) # the samllest is -4.499787

zeros <- which(mod_df$density_SD_log10 == 0)
min(mod_df$density_SD_log10[-zeros]) # the samllest is -4.499787

# so we replace -Inf with -4.5 
# remove zeros
mod_df2 <- mutate(mod_df, abundance_log10 = if_else(abundance_log10 == -Inf, -4.5, abundance_log10),
                  abundance_SD_log10 = if_else(is.na(abundance_SD_log10) == T, 0.51, abundance_SD_log10))

# write a function that will do this for any given species
get_predicted_abundance_function <- function(species) {
  
  species <- enquo(species)
  
  # subset that data to a given species of interest
  dat <- mod_df2 %>%
    dplyr::filter(COMMON_NAME == !!species)
    
  # and filter the detectability model as well
  dat_detect <- detectability %>%
    dplyr::filter(ebird_COMMON_NAME == !!species)
  
  dat2 <- dat %>%
    mutate(predicted=10^(dat_detect$intercept + dat_detect$slope*(dat$abundance_log10)))
  
  return(dat2)
  
}

# test the function with a given species
bcch <- get_predicted_abundance_function("Black-capped Chickadee")

sum(bcch$predicted*bcch$area_square_miles)
sum(bcch$density*bcch$area_square_miles)

ggplot(bcch, aes(x=predicted, y=density))+
  geom_point(aes(size=total_obs))+
  xlab("Predicted density (log10)")+
  ylab("Observed density (log10)")+
  geom_smooth(method="lm", color="orange")+
  scale_x_log10(labels=comma)+
  scale_y_log10(labels=comma)+
  scale_size_continuous(labels=comma)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ggtitle(bcch$COMMON_NAME)+
  labs(size="Total observations")

# do it for three more species (the same species as other example species)
# and then make a four panel figure again
amcr <- get_predicted_abundance_function("American Crow")
amro <- get_predicted_abundance_function("American Robin")
bggn <- get_predicted_abundance_function("Blue-gray Gnatcatcher")

# combine data and then plot
amcr %>%
  bind_rows(amro) %>%
  bind_rows(bcch) %>%
  bind_rows(bggn) %>%
  ggplot(., aes(x=predicted, y=density, color=COMMON_NAME))+
  geom_point(aes(size=total_obs))+
  xlab("Predicted density (log10)")+
  ylab("Observed density (log10)")+
  geom_smooth(method = "lm", formula = y ~ x)+
  facet_wrap( ~ COMMON_NAME, scales = "free")+
  theme_bw()+
  scale_x_log10(labels=comma)+
  scale_y_log10(labels=comma)+
  scale_size_continuous(labels=comma)+
  theme(axis.text=element_text(color="black"))+
  labs(size="Total observations")+
  guides(color=FALSE)+
  scale_color_brewer(palette = "Set2")

ggsave("Figures/Supplementary_figures/example_observed_vs_predicted_density.png", 
       width=6.8, height=6.5, units="in")

# now run the function for every species in the dataset
function_results <- lapply((mod_df2 %>%
                             #dplyr::filter(SEASON != "Wintering") %>%
                             dplyr::select(COMMON_NAME) %>%
                             distinct() %>%
                             .$COMMON_NAME), function(x) {get_predicted_abundance_function(x)})

# put into a df
df_results_from_function <- bind_rows(function_results)

# now summarize to get the estimated abundance
# from the 'observed' values
# to the predicted abundance from our detectability function
summarized_estimates <- df_results_from_function %>%
  group_by(COMMON_NAME) %>%
  summarize(observed_estimate=sum(Population_Estimate),
            predicted_estimate=sum(predicted*area_square_miles))

ggplot(summarized_estimates, aes(x=predicted_estimate, y=observed_estimate))+
  geom_point(color="black")+
  xlab("Predicted population estimate")+
  ylab("Observed population estimate")+
  geom_smooth(method="lm", color="orange")+
  scale_x_log10(labels=comma)+
  scale_y_log10(labels=comma)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_abline(slope=1, intercept=0, color="blue", linetype="dashed", size=2)+
  ggtitle(paste0("N = ", nrow(summarized_estimates), " species"))

ggsave("Figures/Supplementary_figures/observed_vs_predicted_density_of_training_species.png", 
       width=7.8, height=5.9, units="in")


summary(lm(observed_estimate ~ predicted_estimate, data=summarized_estimates))
