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
library(purrr)
library(lme4)
library(MuMIn)

# read in data that was imputed
# necessary for the factors and such within
dat <- readRDS("Intermediate data files/file_for_imputation_15p.RDS")

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
imp <- readRDS("Intermediate data files/imputation_results_15p.RDS")

# try to make a density plot of imputed and observed values
# might not work because 'imp' is really large
# takes forever and doesn't finish on my machine
# Will manually make something similar below
#densityplot(imp)



# get imputed data as a list
imp_list <- mice::complete(imp, "long")

# join factors with imputation results
idat <- imp_list %>% 
  left_join(., factors)

# make plot of imputed vs observed data points
# get 10 imputations to make it more manageable
ten_imps <- idat %>%
  dplyr::filter(.imp %in% sample_n(data.frame(imps=c(1:100)), 10)$imps)

# do this for density and density_se
observed <- ten_imps %>%
  dplyr::filter(Data=="observed") %>%
  distinct(density)

imputed <- ten_imps %>%
  dplyr::filter(Data=="imputed")

density <- ggplot()+
  geom_density(data=imputed, aes(x=density, color="red", group=`.imp`, alpha=0.6), size=0.4)+
  geom_density(data=observed, aes(x=density, color="blue"), linetype="dashed", size=1.2)+
  scale_color_identity(name = "Model fit",
                       breaks = c("red", "blue"),
                       labels = c("imputed", "observed"),
                       guide = "legend")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Density estimate of birds (log10)")+
  ylab("Density of datapoints")+
  guides(linetype=FALSE)+
  guides(alpha=FALSE)+
  guides(size=FALSE)+
  guides(color=FALSE)

density

# do this for density and density_se
observed <- ten_imps %>%
  dplyr::filter(Data=="observed") %>%
  distinct(dens_se)

imputed <- ten_imps %>%
  dplyr::filter(Data=="imputed")

density_se <- ggplot()+
  geom_density(data=imputed, aes(x=dens_se, color="red", group=`.imp`, alpha=0.6), size=0.4)+
  geom_density(data=observed, aes(x=dens_se, color="blue"), linetype="dashed", size=1.2)+
  scale_color_identity(name = "Model fit",
                       breaks = c("red", "blue"),
                       labels = c("imputed", "observed"),
                       guide = "legend")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("SE of density estimate of birds")+
  ylab("")+
  guides(linetype=FALSE)+
  guides(alpha=FALSE)+
  guides(size=FALSE)

density_se

density + density_se + plot_layout(ncol=2)

ggsave("Figures/Supplementary_figures/imputation_check_of_plausible_values_15p.png", width=6.9, height=4.8, units="in")

# Now get the imputed values for each of our held-out species
setwd("Intermediate data files/mice_validation/")
imp_test <- list.files(pattern = ".RDS") %>%
  map(readRDS) %>% 
  bind_rows()
setwd("../..")

# get observed data for the training species
obs_test <- idat %>% 
  dplyr::filter(ebird_COMMON_NAME %in% unique(imp_test$ebird_COMMON_NAME)) %>% 
  dplyr::select(grid_id, ebird_COMMON_NAME, density, dens_se) %>% 
  distinct() %>%
  mutate(Data="observed")

# summarize the imputation data
imp_summary <- imp_test %>%
  ungroup() %>%
  group_by(grid_id, ebird_COMMON_NAME) %>%
  summarize(
    mean_density = mean(density, na.rm = TRUE),
    max_density  = max(density, na.rm = TRUE),
    min_density  = min(density, na.rm = TRUE),
    
    quan_.05_density = quantile(density, 0.05, na.rm = TRUE),
    quan_.95_density = quantile(density, 0.95, na.rm = TRUE),
    
    mean_se = mean(dens_se, na.rm = TRUE),
    max_se  = max(dens_se, na.rm = TRUE),
    min_se  = min(dens_se, na.rm = TRUE),
    
    quan.05_se = quantile(dens_se, 0.05, na.rm = TRUE),
    quan.95_se = quantile(dens_se, 0.95, na.rm = TRUE)
  ) %>%
  ungroup()

testing_data <- imp_summary %>%
  left_join(., obs_test, by=c("grid_id", "ebird_COMMON_NAME"))

# run a lmer to get an R2 for the fit of the
# predicted versus observed estimates
# from the overall leave-one-out mice imputation process
mod_all <- lmer(mean_density ~ density + (1|ebird_COMMON_NAME), data=testing_data)
summary(mod_all)
r.squaredGLMM(mod_all)

# each point represents a grid/species combination
# we have error on the y, but too difficult to plot it here
all_plot <- ggplot(testing_data, aes(x=density, y=mean_density))+
  geom_point()+
  ylab("Mean imputed density (log10)")+
  xlab("Observed density (log10)")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm", color="orange")+
  ggtitle(expression("Species x grid combinations -"~R^2~ "=0.76"))+
  theme(plot.title = element_text(size = 10, face = "bold")) +
  scale_x_continuous(limits = c(-7.5, 5.5)) +
  scale_y_continuous(limits = c(-5.0, 3))

all_plot

# Now look at the same thing
# but collapse one level further
# so that each species has a mean density
# compared to the species mean density
# from the observed data
species_level_summary <- imp_summary %>%
  left_join(., obs_test, by=c("grid_id", "ebird_COMMON_NAME")) %>%
  group_by(ebird_COMMON_NAME) %>%
  summarize(mean_observed_density=mean(density),
            mean_imputed_density=mean(mean_density))

# run a simple lm
# to see the correlation
mod_species <- lm(mean_imputed_density ~ mean_observed_density, data=species_level_summary)
summary(mod_species)


species_plot <- ggplot(species_level_summary, aes(x=mean_observed_density, y=mean_imputed_density))+
  geom_point()+
  ylab("")+
  xlab("Observed density (log10)")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm", color="orange")+
  ggtitle(expression("Mean species across grids -"~R^2~ "=0.59"))+
  theme(plot.title = element_text(size = 10, face = "bold")) +
  scale_x_continuous(limits = c(-3, 2.5)) +
  scale_y_continuous(limits = c(-3.5, 1))

species_plot


all_plot + species_plot + plot_layout(ncol=2)

ggsave("Figures/Supplementary_figures/mice_imputation_validation_plot_15p.png", width=6.9, height=4.8, units="in")


# some quick quantifying of stuff
# first the percentage of observations that
# were within the range of imputated data
testing_data %>%
  mutate(in_range=ifelse(density>=min_density & density<=max_density, "Yes", "No")) %>%
  group_by(in_range) %>%
  summarize(N=n()) %>%
  mutate(total=nrow(testing_data)) %>%
  mutate(percent=(N/total)*100)

# now the percent of observations that were within the
# 95% range of data
testing_data %>%
  mutate(in_range=ifelse(density>=quan_.05_density & density<=quan_.95_density, "Yes", "No")) %>%
  group_by(in_range) %>%
  summarize(N=n()) %>%
  mutate(total=nrow(testing_data)) %>%
  mutate(percent=(N/total)*100)

# another way of visualizing
# for just 12 of the species
testing_data %>%
  dplyr::filter(ebird_COMMON_NAME %in% c(testing_data %>%
                                           dplyr::select(ebird_COMMON_NAME) %>%
                                           distinct() %>%
                                           sample_n(12) %>%
                                           .$ebird_COMMON_NAME)) %>%
  ggplot(.)+
  geom_linerange(aes(x=as.factor(grid_id), ymin=min_density, ymax=max_density), size=0.6)+
  facet_wrap(~ebird_COMMON_NAME, scales="free", ncol=3)+
  coord_flip()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_point(aes(x=as.factor(grid_id), y=density), color="red")+
  ylab("Density of birds (log10)")+
  xlab("grid_id")+
  theme(strip.text.x = element_text(size = 8))+
  theme(axis.text.y=element_blank())

ggsave("Figures/Supplementary_figures/mice_imputation_species_example_15p.png", width=7.4, height=8.6)

# calculate the number of grids a training species is in
grids_per_species <- testing_data %>%
  group_by(ebird_COMMON_NAME) %>%
  summarize(N=n())

min(grids_per_species$N)
max(grids_per_species$N)
mean(grids_per_species$N)
sd(grids_per_species$N)

grids_training <- ggplot(grids_per_species, aes(x=N))+
  geom_histogram(color="black", fill="lightblue", bins=50)+
  xlab("Number of grids a species occupies")+
  ylab("Number of species")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ggtitle("A")

grids_training

ggsave("Figures/Supplementary_figures/grids_per_species_training_15p.png", width=5.6, height=4.9, units="in")


# calculate the number of grids all species are in
grids_per_species_all <- factors %>%
  group_by(ebird_COMMON_NAME) %>%
  summarize(N=n())

min(grids_per_species_all$N)
max(grids_per_species_all$N)
mean(grids_per_species_all$N)
sd(grids_per_species_all$N)

grids_all <- ggplot(grids_per_species_all, aes(x=N))+
  geom_histogram(color="black", fill="lightblue", bins=50)+
  xlab("Number of grids a species occupies")+
  ylab("Number of species")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  ggtitle("B")

grids_all

ggsave("Figures/Supplementary_figures/grids_per_species_all_15p.png", width=5.6, height=4.9, units="in")

grids_training + grids_all + plot_layout(ncol=1)


ggsave("Figures/Supplementary_figures/grids_per_species_15p.png", width=5.6, height=7.9, units="in")
