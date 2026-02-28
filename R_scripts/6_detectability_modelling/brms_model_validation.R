# This is an R script to test different ways of including data in the original
# detectability model (brms model, which was previously MCMCglmm model)

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
               forcats
)


# running rstan + brms
options(mc.cores = 24)
rstan_options(auto_write = TRUE)

# read in the data for detectability modelling
mod_df <- readRDS("Intermediate data files/data_for_detectability_modelling.RDS")

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
  theme_classic()+
  geom_smooth(method="lm", color="orange")+
  theme(axis.text=element_text(color="black"))+
  xlab("eBird relative abundance (log10)")+
  ylab("Observed density (log10)")

# First read in the brms_mod that is the big one with 10000 iterations and 4 chains etc
# and do some summary stuff built in from brms and loo
# we will want to show that we don't lose a lot of information
# by repeating the model with a much smaller number of iterations
# just as a sanity check
brms_mod <- readRDS(file = "./Intermediate data files/brms_mod.rds")

big_mod_loo <- loo(brms_mod)

big_mod_loo

plot(big_mod_loo)

bayesian_R2_mean <- bayes_R2(brms_mod)
bayesian_R2_mean
bayesian_R2_median <- bayes_R2(brms_mod, robust=TRUE)

loosian_R2 <- loo_R2(brms_mod)
loosian_R2

# Get the species-specific summary stuff here
# Summarize data for this model
means <- fixef(brms_mod)
blups<-ranef(brms_mod)

intercepts <- means[1,1] + as.numeric(blups$species[ , 1, 1])
slopes <-   means[2,1] + + as.numeric(blups$species[ , 1, 2])

# errors
intercepts_se <- sqrt(means[1,2]^2 + as.numeric(blups$species[ , 2, 1])^2)
slopes_se <-   sqrt(means[2,2]^2 + as.numeric(blups$species[ , 2, 2])^2)

n_spp <- length(unique(mod_df2$species))
blups_df_big <- tibble(species = rep(attr(blups$species, "dimnames")[[1]], 2),
                       type = rep(c("intercept","slope"), each = n_spp),
                       estimate = c(intercepts, slopes), 
                       se =  c(intercepts_se, slopes_se))

# replot the data using the brms mean intercept and slope
ggplot(mod_df2, aes(y=density_log10, x=abundance_log10))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("eBird relative abundance (log10)")+
  ylab("Observed density (log10)")+
  geom_abline(intercept = fixef(brms_mod)[1], 
              slope     = fixef(brms_mod)[2],
              color="orange", linetype="solid", size=2.3)

post_final <- posterior_samples(brms_mod)

ggplot(mod_df2, aes(y=density_log10, x=abundance_log10))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("eBird relative abundance (log10)")+
  ylab("Observed density (log10)")+
  geom_abline(intercept = post_final[1:20, 1], 
            slope     = post_final[1:20, 2],
            size = 1.4, color="orange", alpha = .3)


# this is our main brms model used for detectability
# but although we use 10000 iterations in the
# big one and four chains
# here I'll use 2 chains and 4000 iterations to do comparisons with
# different things
brms_mod_final <- brm(density_log10 ~ 1 + me(abundance_log10, abundance_SD_log10) + 
                  (1 + me(abundance_log10, abundance_SD_log10) |species), 
                data = mod_df2, save_mevars = T, chains = 2, iter = 4000, warmup = 2000)


# Summarize data for this model
means <- fixef(brms_mod_final)
blups<-ranef(brms_mod_final)

intercepts <- means[1,1] + as.numeric(blups$species[ , 1, 1])
slopes <-   means[2,1] + + as.numeric(blups$species[ , 1, 2])

# errors
intercepts_se <- sqrt(means[1,2]^2 + as.numeric(blups$species[ , 2, 1])^2)
slopes_se <-   sqrt(means[2,2]^2 + as.numeric(blups$species[ , 2, 2])^2)

n_spp <- length(unique(mod_df2$species))
blups_df_final <- tibble(species = rep(attr(blups$species, "dimnames")[[1]], 2),
                   type = rep(c("intercept","slope"), each = n_spp),
                   estimate = c(intercepts, slopes), 
                   se =  c(intercepts_se, slopes_se))

# Is the "BIG" model similar to one with a smaller number of
# iterations
estimate_check <- blups_df_big %>%
  left_join(blups_df_final, by=c("species", "type")) %>%
  rename(estimate_big=estimate.x,
         se_big=se.x,
         estimate_small=estimate.y,
         se_small=se.y) %>%
  ggplot(., aes(x=estimate_big, y=estimate_small))+
  geom_point()+
  geom_smooth()+
  theme_bw()+
  facet_wrap(~type, scales="free")+
  xlab("Estimate from model with 10000 iterations and 4 chains")+
  ylab("Estimate from model with 4000 iterations and 2 chains")

estimate_check

se_check <- blups_df_big %>%
  left_join(blups_df_final, by=c("species", "type")) %>%
  rename(estimate_big=estimate.x,
         se_big=se.x,
         estimate_small=estimate.y,
         se_small=se.y) %>%
  ggplot(., aes(x=se_big, y=se_small))+
  geom_point()+
  geom_smooth()+
  theme_bw()+
  facet_wrap(~type, scales="free")+
  xlab("SE from model with 10000 iterations and 4 chains")+
  ylab("SE from model with 4000 iterations and 2 chains")

se_check

estimate_check + se_check + plot_layout(ncol=1)

# Okay, so this basically validates the idea of testing different things with brms
# on a smaller model with smaller chains and iterations
# means can more-or-less ignore the warnings about ESS etc.

# Now run a brms model without the zeros
mod_df3 <- mod_df %>%
  dplyr::filter(abundance_log10 != -Inf) %>%
  dplyr::filter(complete.cases(abundance_SD_log10))

sum(is.na(mod_df3$abundance_SD_log10))
# checking both of these (abundance_log10 and abundance_SD_log10) how it looks

hist(mod_df3$abundance_log10)
hist(log(mod_df3$abundance_SD_log10))

ggplot(mod_df3, aes(y=density_log10, x=abundance_log10))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm", color="orange")+
  theme(axis.text=element_text(color="black"))+
  xlab("eBird relative abundance (log10)")+
  ylab("Observed density (log10)")


# this is our main brms model used for detectability
# but although we use 10000 iterations in the
# big one and four chains
# here I'll use 2 chains and 4000 iterations to do comparisons with
# different things
# now run a model without the zeros as a reviewer mentioned this
brms_mod_no_abund_zeros <- brm(density_log10 ~ 1 + me(abundance_log10, abundance_SD_log10) + 
                        (1 + me(abundance_log10, abundance_SD_log10) |species), 
                      data = mod_df3, save_mevars = T, chains = 2, iter = 4000, warmup = 2000)


# Summarize data for this model
means <- fixef(brms_mod_no_abund_zeros)
blups<-ranef(brms_mod_no_abund_zeros)

intercepts <- means[1,1] + as.numeric(blups$species[ , 1, 1])
slopes <-   means[2,1] + + as.numeric(blups$species[ , 1, 2])

# errors
intercepts_se <- sqrt(means[1,2]^2 + as.numeric(blups$species[ , 2, 1])^2)
slopes_se <-   sqrt(means[2,2]^2 + as.numeric(blups$species[ , 2, 2])^2)

n_spp <- length(unique(mod_df3$species))
blups_df_no_abund_zeros <- tibble(species = rep(attr(blups$species, "dimnames")[[1]], 2),
                         type = rep(c("intercept","slope"), each = n_spp),
                         estimate = c(intercepts, slopes), 
                         se =  c(intercepts_se, slopes_se))

# replot the data using the brms mean intercept and slope
ggplot(mod_df3, aes(y=density_log10, x=abundance_log10))+
  geom_point()+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  xlab("eBird relative abundance (log10)")+
  ylab("Observed density (log10)")+
  geom_abline(intercept = fixef(brms_mod_no_abund_zeros)[1], 
              slope     = fixef(brms_mod_no_abund_zeros)[2],
              color="orange", linetype="solid", size=2.3)

post_no_abund_zeros <- posterior_samples(brms_mod_no_abund_zeros)

ggplot(mod_df3, aes(y=density_log10, x=abundance_log10))+
  geom_point()+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  xlab("eBird relative abundance (log10)")+
  ylab("Observed density (log10)")+
  geom_abline(intercept = post_no_abund_zeros[1:20, 1], 
              slope     = post_no_abund_zeros[1:20, 2],
              size = 1.4, color="orange", alpha = .3)

# See if by including the zeros in the model
# whether or not they are influencing
# the fits respective for each species etc.
# this is only for a smaller subset of species
# so won't be the same number of species for each
effect_of_zeros <- blups_df_final %>%
  right_join(blups_df_no_abund_zeros, by=c("species", "type")) %>%
  rename(estimate_final=estimate.x,
         se_final=se.x,
         estimate_no_zeros=estimate.y,
         se_no_zeros=se.y)

estimates <- ggplot(effect_of_zeros, aes(x=estimate_final, y=estimate_no_zeros))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()+
  facet_wrap(~type, scales="free")+
  xlab("Estimate with all data")+
  ylab("Estimate from model with no zeros for eBird")+
  theme(axis.title=element_text(size=10))+
  ggtitle(paste0("N=", length(unique(effect_of_zeros$species)), " species"))

estimates

ses <- ggplot(effect_of_zeros, aes(x=se_final, y=se_no_zeros))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()+
  facet_wrap(~type, scales="free")+
  xlab("SE with all data")+
  ylab("SE from model with no zeros for eBird")+
  theme(axis.title=element_text(size=10))

ses

estimates + ses + plot_layout(ncol=1)

effect_of_zeros %>%
  dplyr::filter(type=="intercept") %>%
  lm(estimate_no_zeros ~ estimate_final, data=.) %>%
  summary()

effect_of_zeros %>%
  dplyr::filter(type=="slope") %>%
  lm(estimate_no_zeros ~ estimate_final, data=.) %>%
  summary()

effect_of_zeros %>%
  dplyr::filter(type=="intercept") %>%
  lm(se_no_zeros ~ se_final, data=.) %>%
  summary()

effect_of_zeros %>%
  dplyr::filter(type=="slope") %>%
  lm(se_no_zeros ~ se_final, data=.) %>%
  summary()

# NOw do another test incorporating error on the y (outcome)
# how many 0 in density_SD_log10
length(which(mod_df$density_SD_log10 == 0))
length(which(is.na(mod_df$density_SD_log10) == T))

# how many NA in density_SD_log10

mod_df4 <- mod_df2 %>% 
  filter(density_SD_log10 != 0 | is.na(density_SD_log10) == F) %>%
  mutate(density_SD_log10=ifelse(density_SD_log10==0, 0.0001, density_SD_log10))


brms_mod_error_on_y <- brm(density_log10 | se(density_SD_log10, sigma = TRUE) ~ 1 + me(abundance_log10, abundance_SD_log10) + 
                                 (1 + me(abundance_log10, abundance_SD_log10) |species), 
                               data = mod_df4, save_mevars = T, chains = 2, iter = 4000, warmup = 2000)

# Summarize data for this model
means <- fixef(brms_mod_error_on_y)
blups<-ranef(brms_mod_error_on_y)

intercepts <- means[1,1] + as.numeric(blups$species[ , 1, 1])
slopes <-   means[2,1] + + as.numeric(blups$species[ , 1, 2])

# errors
intercepts_se <- sqrt(means[1,2]^2 + as.numeric(blups$species[ , 2, 1])^2)
slopes_se <-   sqrt(means[2,2]^2 + as.numeric(blups$species[ , 2, 2])^2)

n_spp <- length(unique(mod_df4$species))
blups_df_error_on_y <- tibble(species = rep(attr(blups$species, "dimnames")[[1]], 2),
                                  type = rep(c("intercept","slope"), each = n_spp),
                                  estimate = c(intercepts, slopes), 
                                  se =  c(intercepts_se, slopes_se))

# replot the data using the brms mean intercept and slope
ggplot(mod_df4, aes(y=density_log10, x=abundance_log10))+
  geom_point()+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  xlab("eBird relative abundance (log10)")+
  ylab("Observed density (log10)")+
  geom_abline(intercept = fixef(brms_mod_error_on_y)[1], 
              slope     = fixef(brms_mod_error_on_y)[2],
              color="orange", linetype="solid", size=2.3)

post_error_on_y <- posterior_samples(brms_mod_error_on_y)

ggplot(mod_df4, aes(y=density_log10, x=abundance_log10))+
  geom_point()+
  theme_classic()+
  theme(axis.text=element_text(color="black"))+
  xlab("eBird relative abundance (log10)")+
  ylab("Observed density (log10)")+
  geom_abline(intercept = post_error_on_y[1:20, 1], 
              slope     = post_error_on_y[1:20, 2],
              size = 1.4, color="orange", alpha = .3)

# check that the intercept and slope are not biased
# See if by including the zeros in the model
# whether or not they are influencing
# the fits respective for each species etc.
# this is only for a smaller subset of species
# so won't be the same number of species for each
effect_of_measurement_error_on_y <- blups_df_final %>%
  right_join(blups_df_error_on_y, by=c("species", "type")) %>%
  rename(estimate_final=estimate.x,
         se_final=se.x,
         estimate_measurement_error_on_y=estimate.y,
         se_measurement_error_on_y=se.y)

estimates <- ggplot(effect_of_measurement_error_on_y, aes(x=estimate_final, 
                                                          y=estimate_measurement_error_on_y))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()+
  facet_wrap(~type, scales="free")+
  xlab("Estimate with all data")+
  ylab("Estimate from model with measurement error on y")+
  theme(axis.title=element_text(size=10))+
  ggtitle(paste0("N=", length(unique(effect_of_measurement_error_on_y$species)), " species"))

estimates

effect_of_measurement_error_on_y %>%
  dplyr::filter(type=="intercept") %>%
  lm(estimate_measurement_error_on_y ~ estimate_final, data=.) %>%
  summary()

effect_of_measurement_error_on_y %>%
  dplyr::filter(type=="slope") %>%
  lm(estimate_measurement_error_on_y ~ estimate_final, data=.) %>%
  summary()

# compare the standard error for
# a model with all species and error on the x only
# and a model with error on the x and y
error_comparison <- blups_df_final %>%
  right_join(., blups_df_error_on_y, by=c("species", "type")) %>%
  rename(estimate_full_model=estimate.x,
         se_error_full_model=se.x,
         estimate_measurement_error_on_y=estimate.y,
         se_measurement_error_on_y=se.y) 

slope_error <- error_comparison %>%
  dplyr::filter(type=="slope") %>%
  arrange(se_error_full_model) %>%
  dplyr::select(species, se_error_full_model, se_measurement_error_on_y) %>%
  pivot_longer(!species, names_to="error", values_to="value") %>%
  arrange(value)


slope_error1 <- ggplot(slope_error, aes(x=fct_inorder(species), y=value, color=error))+
  geom_point()+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.y=element_blank())+
  theme(panel.grid.major=element_blank())+
  theme(panel.grid.minor=element_blank())+
  theme(axis.text=element_text(color="black"))+
  ylab("SE of slope")+
  xlab("Species")+
  theme(axis.text.y=element_blank())+
  scale_color_brewer(palette = "Set2")+
  ggtitle(paste0("N=", length(unique(slope_error$species)), " species"))

slope_error1

slope_error2 <- ggplot(slope_error, aes(x=value, fill=error))+
  geom_density(alpha=0.5)+
  theme_bw()+
  xlab("SE of slope")+
  ylab("Density of observations")+
  scale_fill_brewer(palette = "Set2")+
  ggtitle(paste0("N=", length(unique(slope_error$species)), " species"))

slope_error2

intercept_error <- error_comparison %>%
  dplyr::filter(type=="intercept") %>%
  arrange(se_error_full_model) %>%
  dplyr::select(species, se_error_full_model, se_measurement_error_on_y) %>%
  pivot_longer(!species, names_to="error", values_to="value") %>%
  arrange(value)


intercept_error1 <- ggplot(intercept_error, aes(x=fct_inorder(species), y=value, color=error))+
  geom_point()+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.y=element_blank())+
  theme(panel.grid.major=element_blank())+
  theme(panel.grid.minor=element_blank())+
  theme(axis.text=element_text(color="black"))+
  ylab("SE of intercept")+
  xlab("Species")+
  theme(axis.text.y=element_blank())+
  scale_color_brewer(palette = "Set2")

intercept_error1

intercept_error2 <- ggplot(slope_error, aes(x=value, fill=error))+
  geom_density(alpha=0.5)+
  theme_bw()+
  xlab("SE of intercept")+
  ylab("Density of observations")+
  scale_fill_brewer(palette = "Set2")

intercept_error2


# two ways to look at this
slope_error1 + intercept_error1 + plot_layout(ncol=1)

slope_error2 + intercept_error2 + plot_layout(ncol=1)


# Make a plot showing the fixed effects for all three approaches
# and showing that they are largely similar
# we aren't particularly interested in the fixed effects per se
# but as the Reviewer points out, it is and by looking at the fixed effects
# we can show that the model is robust to the different decisions regarding data
# to be included in the final model
ggplot()+
  geom_point(data=mod_df4, aes(y=density_log10, x=abundance_log10), color="transparent")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("eBird relative abundance (log10)")+
  ylab("Observed density (log10)")+
  geom_abline(aes(intercept = post_error_on_y[1:20, 1], 
              slope     = post_error_on_y[1:20, 2],
              color="orange"), size = 1.4, alpha = .3)+
  geom_abline(aes(intercept = post_no_abund_zeros[1:20, 1], 
              slope     = post_no_abund_zeros[1:20, 2],
              color="blue"), size = 1.4, alpha = .3)+
  geom_abline(aes(intercept = post_final[1:20, 1], 
              slope     = post_final[1:20, 2],
              color="green"), size = 1.4, alpha = .3)+
  scale_color_identity(name = "Model:",
                       breaks = c("orange", "blue", "green"),
                       labels = c("measurement_error_on_y", "no_zeros_included", "full_model"),
                       guide = "legend")+
  theme(legend.position="bottom")





