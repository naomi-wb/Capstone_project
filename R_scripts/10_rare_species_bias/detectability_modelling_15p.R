# this script is used to fit a model
# which estimates the 'citizen science detectability'
# for all the 'training data' we have data for
# the detectability of a species is then a function of
# the slope and intercept of this model fit



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
               tidybayes
)
# running rstan + brms
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# read in the data for detectability modelling
mod_df <- readRDS("Intermediate data files/data_for_detectability_modelling_15p.RDS")

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
  
# checking both of these (abundance_log10 and abundance_SD_log10) how it looks

hist(mod_df2$abundance_log10)
hist(log(mod_df2$abundance_SD_log10))

# plot the relationship between density and abundance
# for all data points in our dataset
ggplot(mod_df2, aes(y=density_log10, x=abundance_log10))+
  geom_point()+
  theme_bw()+
  geom_smooth(method="lm", color="orange")+
  theme(axis.text=element_text(color="black"))+
  xlab("eBird relative abundance (log10)")+
  ylab("Observed density (log10)")

# make a figure of the number of data points per species in
# the detectability model
mod_df2 %>%
  group_by(COMMON_NAME) %>%
  summarize(N=n()) %>%
  ggplot(., aes(x=N))+
  geom_histogram(bins=50, color="black", fill="lightblue")+
  scale_x_log10(labels=comma)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Number of unique data points")+
  ylab("Number of species")

ggsave("Figures/Supplementary_figures/number_of_data_points_per_species_15p.png",
       width=6.8, height=5.6, units="in")


# brms model with measurement error on abundance_log10!!
# I guess if we want we could get error for density_log10 as well
brms_mod <- brm(density_log10 ~ 1 + me(abundance_log10, abundance_SD_log10) + 
                  (1 + me(abundance_log10, abundance_SD_log10) |species), 
                data = mod_df2, save_mevars = T,chains = 4, iter = 10000, warmup = 2000)

#saveRDS(brms_mod,  file = "./Intermediate data files/brms_mod.rds")
#brms_mod <- readRDS(file = "./Intermediate data files/brms_mod.rds")
brms_mod_15p <- readRDS(file = "./Intermediate data files/brms_mod_15p.rds")

summary(brms_mod_15p)
means <- fixef(brms_mod_15p)
blups<-ranef(brms_mod_15p)

intercepts <- means[1,1] + as.numeric(blups$species[ , 1, 1])
slopes <-   means[2,1] + + as.numeric(blups$species[ , 1, 2])

# errors
intercepts_se <- sqrt(means[1,2]^2 + as.numeric(blups$species[ , 2, 1])^2)
slopes_se <-   sqrt(means[2,2]^2 + as.numeric(blups$species[ , 2, 2])^2)

n_spp <- length(unique(mod_df2$species))
blups_df <- tibble(species = rep(attr(blups$species, "dimnames")[[1]], 2),
                        type = rep(c("intercept","slope"), each = n_spp),
                        estimate = c(intercepts, slopes), 
                        se =  c(intercepts_se, slopes_se)) #

#hist(blups_df$estimate[blups_df$type == "intercept"],breaks = 20) 
#hist(blups_df$estimate[blups_df$type == "slope"],breaks = 20) 
#hist(blups_df$se[blups_df$type == "intercept"],breaks = 20) 
#hist(blups_df$se[blups_df$type == "slope"],breaks = 20) 

a <- blups_df %>%
  dplyr::filter(type=="intercept") %>%
  ggplot(., aes(x=estimate))+
  geom_histogram(bins=50, color="black", fill="lightblue")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("brms intercepts")+
  ylab("Number of species")

a

b <- blups_df %>%
  dplyr::filter(type=="intercept") %>%
  ggplot(., aes(x=se))+
  geom_histogram(bins=50, color="black", fill="lightblue")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("brms intercept se")+
  ylab("Number of species")

b

c <- blups_df %>%
  dplyr::filter(type=="slope") %>%
  ggplot(., aes(x=estimate))+
  geom_histogram(bins=50, color="black", fill="lightblue")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("brms slopes")+
  ylab("Number of species")

c

d <- blups_df %>%
  dplyr::filter(type=="slope") %>%
  ggplot(., aes(x=se))+
  geom_histogram(bins=50, color="black", fill="lightblue")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("brms slopes se")+
  ylab("Number of species")

d

a + b + c + d + plot_layout(ncol=2)

ggsave("Figures/Supplementary_figures/brms_histograms_15p.png",
       width=6.8, height=5.6, units="in")

spp_slope <- blups_df %>% filter(type == "slope") %>% select(species, estimate, se)
arrange(spp_slope, estimate)
write_csv(spp_slope, "./Intermediate data files/spp_slope_15p.csv")
write_csv(blups_df, "./Intermediate data files/spp_all_15p.csv")

# get posterior draws (N=20) and plot on the raw data points
# for a supplementary figure
post_final <- posterior_samples(brms_mod_15p)

ggplot(mod_df2, aes(y=density_log10, x=abundance_log10))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("eBird relative abundance (log10)")+
  ylab("Observed density (log10)")+
  geom_abline(intercept = post_final[1:20, 1], 
              slope     = post_final[1:20, 2],
              size = 1.4, color="orange", alpha = .3)

ggsave("Figures/Supplementary_figures/training_data_density_vs_abundance_15p.png",
       width=6.8, height=5.6, units="in")
