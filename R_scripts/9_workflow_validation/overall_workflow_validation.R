# An assessment of overall workflow validation
# this first relies on the results from
# the script: "brms_and_mice_validation_combined.R"
# which randomly witheld a species and then repeated
# entire analysis without that species
# but wrote out the results for THAT species'
# population abundance distribution when that species
# was not considered in the analysis
# here we get this list of species
# and get their abundance distributions from the analysis
# when they WERE included (when all species were included)
# and compare to see if they are strongly correlated
# as training data and witheld data

# packages
library(readr)
library(dplyr)
library(purrr)
library(HDInterval)
library(scales)
library(ggplot2)

# first set up the workflow from "summarize_imputations2.R"
# which will reproduce the workflow for when all species
# were included
# read in clean clements
# which will come in handy above
clements <- read_csv("Data/clements_clean.csv")

# read in ecological niche and grouping data
niches <- read_csv("Data/ecological_niche_assignment/41559_2019_1070_MOESM3_ESM.csv") %>%
  dplyr::rename(TipLabel=Binomial)

# read in data that was imputed
# necessary for the factors and such within
dat <- readRDS("Intermediate data files/file_for_imputation.RDS")

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
imp <- readRDS("Intermediate data files/imputation_results.RDS")

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

# want to write out the data to make the distributions for a 
# select number of species
# which can be presented in the main manuscript
# save each species
# as a separate RDS and then they can be joined later on
example_species_function <- function(example_species) {
  
  species_col_number <- sppdat %>%
    mutate(row_id=1:nrow(.)) %>%
    dplyr::filter(ebird_COMMON_NAME == example_species) %>%
    .$row_id
  
  species_dat <- tibble(non = spp_num_non[[species_col_number]],
                        imp = spp_num_imp[[species_col_number]]) %>%
    mutate(ebird_COMMON_NAME=sppdat$ebird_COMMON_NAME[[species_col_number]])
  
  return(species_dat)
  
}

# now read in the example species from the list of species 
# which were withheld from the analysis process

example_species <- gsub(".RDS", "", list.files("Intermediate data files/workflow_validation/with_held/"))


full_data_results <- bind_rows(lapply(example_species, example_species_function)) %>%
  mutate(Data="All species")


# Now read in the results from when a species was withheld from the analysis
setwd("Intermediate data files/workflow_validation/with_held/")
with_held_data_results <- list.files(pattern = ".RDS") %>%
  map(readRDS) %>% 
  bind_rows() %>%
  mutate(Data="Withheld")

setwd("../../..")


# plot the median population estimate versus the median population estimate
full_summary <- full_data_results %>%
  group_by(Data, ebird_COMMON_NAME) %>%
  summarize(non_median=median(non),
            non_.05_quan=quantile(non, .05),
            non_.95_quan=quantile(non, .95),
            imp_median=median(imp),
            imp_.05_quan=quantile(imp, .05),
            imp_.95_quan=quantile(imp, .95))

with_held_summary <- with_held_data_results %>%
  group_by(Data, ebird_COMMON_NAME) %>%
  summarize(non_median=median(non),
            non_.05_quan=quantile(non, .05),
            non_.95_quan=quantile(non, .95),
            imp_median=median(imp),
            imp_.05_quan=quantile(imp, .05),
            imp_.95_quan=quantile(imp, .95))

# without_imp_error <- full_summary %>%
#   left_join(., with_held_summary, by="ebird_COMMON_NAME") %>%
#   ggplot(., aes(x=non_median.x, y=non_median.y))+
#   geom_point()+
#   xlab("Abundance from full model")+
#   ylab("Abundance when a species is withheld")+
#   scale_x_log10(labels=comma)+
#   scale_y_log10(labels=comma)+
#   theme_bw()+
#   theme(axis.text=element_text(color="black"))+
#   geom_smooth(method="lm", color="orange")+
#   ggtitle("Without imputation error")
# 
# without_imp_error

with_imp_error <- full_summary %>%
  left_join(., with_held_summary, by="ebird_COMMON_NAME") %>%
  ggplot(., aes(x=imp_median.x, y=imp_median.y))+
  geom_point()+
  xlab("Abundance from full model")+
  ylab("Abundance when species is withheld")+
  scale_x_log10(labels=comma)+
  scale_y_log10(labels=comma)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm", color="orange")

with_imp_error

# get R2 of this relationship
full_summary %>%
  left_join(., with_held_summary, by="ebird_COMMON_NAME") %>%
  lm(log10(imp_median.x) ~ log10(imp_median.y), data=.) %>%
  summary()

# remake plot with R2 of this relationship!
full_summary %>%
  left_join(., with_held_summary, by="ebird_COMMON_NAME") %>%
  ggplot(., aes(x=imp_median.x, y=imp_median.y))+
  geom_point()+
  xlab("Abundance from full model")+
  ylab("Abundance when a species is withheld")+
  scale_x_log10(labels=comma)+
  scale_y_log10(labels=comma)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm", color="orange")+
  ggtitle(expression("N=684 species -"~R^2~ "=0.94"))

ggsave("Figures/Supplementary_figures/overall_workflow_validation.png", width=5.6, height=4.8, units="in")

# pick random 10 species
# and show geom_density
# differences between imputed and training data
ten_species <- example_species %>%
  as.data.frame() %>%
  sample_n(10) %>%
  .[[1]]

full_data_results %>%
  bind_rows(with_held_data_results) %>%
  dplyr::filter(ebird_COMMON_NAME %in% ten_species) %>%
  ggplot(., aes(x=non, fill=Data))+
  geom_density(alpha=0.5)+
  scale_fill_brewer(palette = "Set1")+
  theme_bw()+
  scale_x_log10(labels=comma)+
  theme(axis.text=element_text(color="black"))+
  facet_wrap(~ebird_COMMON_NAME, scales="free_x", ncol=2)+
  xlab("Abundance (log10)")+
  ggtitle("Without imputation error")

full_data_results %>%
  bind_rows(with_held_data_results) %>%
  dplyr::filter(ebird_COMMON_NAME %in% ten_species) %>%
  mutate(Data=case_when(Data=="All species" ~ "Full model",
                        Data=="Withheld" ~ "Species withheld")) %>%
  ggplot(., aes(x=imp, fill=Data))+
  geom_density(alpha=0.5)+
  scale_fill_brewer(palette = "Set1")+
  theme_bw()+
  scale_x_log10(labels=comma)+
  theme(axis.text=element_text(color="black"))+
  facet_wrap(~ebird_COMMON_NAME, scales="free_x", ncol=2)+
  xlab("Number of individual birds")+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())+
  theme(axis.text.x=element_text(hjust=0.7))+
  theme(axis.ticks.x=element_blank())
  
ggsave("Figures/Supplementary_figures/overall_workflow_validation_example_species.png",
       width=8.8, height=8.9, units="in")
