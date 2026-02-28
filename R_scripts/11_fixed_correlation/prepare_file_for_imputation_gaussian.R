####################################
# prepare_file_for_imputation (REFACTORED + FIXED)
####################################

####################################
# 0. PACKAGE SETUP
####################################

cran_packages <- c(
  "tidyverse", "magrittr", "broom", "lme4", "GGally",
  "corrplot", "here", "rapportools", "mvtnorm", "sf",
  "ggcorrplot", "posterior", "brms", "rstan"
)

installed <- rownames(installed.packages())
to_install <- setdiff(cran_packages, installed)
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(mvtnorm)
  library(sf)
  library(brms)
  library(posterior)
})

####################################
# 1. READ STATIC DATA
####################################

clements_clean <- read_csv("Data/clements_clean.csv")

body_size <- read_csv("Data/body size data/cleaned_body_size_data.csv") %>%
  select(3, 7:9)

color <- read_csv("Data/color data/color_data_by_species.csv")

IUCN <- read_csv("Data/IUCN categories/cleaned_IUCN_data.csv")

flock_size <- readRDS("Data/Flock size/flock_size_per_month_for_each_species.RDS") %>%
  group_by(COMMON_NAME) %>%
  summarise(mean_max_flock_size = mean(max_abund), .groups = "drop") %>%
  rename(ebird_COMMON_NAME = COMMON_NAME)

detectability <- read_csv("Intermediate data files/spp_all.csv") %>%
  pivot_wider(names_from = type, values_from = c(estimate, se)) %>%
  rename(
    ebird_COMMON_NAME = species,
    intercept = estimate_intercept,
    slope = estimate_slope,
    intercept_se = se_intercept,
    slope_se = se_slope
  )

####################################
# 2. PREPARE IMPUTATION FILE
####################################

imputation_file <- clements_clean %>%
  distinct(ebird_COMMON_NAME, .keep_all = TRUE) %>%
  select(ebird_COMMON_NAME, TipLabel) %>%
  left_join(flock_size, by = "ebird_COMMON_NAME") %>%
  left_join(color, by = "TipLabel") %>%
  left_join(body_size, by = "TipLabel") %>%
  left_join(IUCN, by = "TipLabel") %>%
  group_by(ebird_COMMON_NAME) %>%
  slice(1) %>%
  mutate(
    IUCN_ordinal = case_when(
      IUCN_category %in% c("EW","EX","CR (PE)","CR (PEW)") ~ 0,
      IUCN_category == "CR" ~ 1,
      IUCN_category == "EN" ~ 2,
      IUCN_category == "VU" ~ 3,
      IUCN_category == "NT" ~ 4,
      IUCN_category == "LC" ~ 5
    ),
    adult_body_mass_g = ifelse(adult_body_mass_g == -999, NA, adult_body_mass_g),
    mass_log10 = log10(adult_body_mass_g),
    flock_size_log10 = log10(mean_max_flock_size),
    species = ebird_COMMON_NAME
  ) %>%
  left_join(detectability, by = "ebird_COMMON_NAME") %>%
  left_join(
    clements_clean %>%
      select(ebird_COMMON_NAME, ebird_SCIENTIFIC_NAME, order, family) %>%
      distinct(),
    by = "ebird_COMMON_NAME"
  ) %>%
  ungroup()

####################################
# 3. READ AND COMBINE EBIRD GRID FILES
####################################

months <- tolower(month.name)

combined_dat <- bind_rows(
  lapply(months, function(m) {
    readRDS(
      paste0(
        "Intermediate data files/abundance_grid_results/",
        m, ".RDS"
      )
    )
  })
)

####################################
# 4. SPECIES × GRID TEMPORAL SUMMARY
####################################

dat <- combined_dat %>%
  filter(mean_abund > 0) %>%
  group_by(grid_id, ebird_COMMON_NAME) %>%
  summarise(
    mean_abund = mean(mean_abund),
    number_months = n_distinct(MONTH),
    mean_species_checklists = mean(obs_with_abund),
    total_species_checklists = sum(obs_with_abund),
    .groups = "drop"
  ) %>%
  mutate(log10_mean_abund = log10(mean_abund))

####################################
# 5. GRID-LEVEL SUMMARY  
####################################

grid_summary <- combined_dat %>%
  filter(mean_abund > 0) %>%
  select(grid_id, total_checklists, area_square_miles) %>%
  distinct() %>%
  group_by(grid_id) %>%
  summarise(
    mean_grid_checklists = mean(total_checklists),
    total_grid_checklists = sum(total_checklists),
    area_square_miles = first(area_square_miles),
    .groups = "drop"
  )

dat <- dat %>%
  left_join(grid_summary, by = "grid_id")

rm(combined_dat)

####################################
# 6. JOIN TRAITS & COMPUTE DENSITY
####################################

dat_imp <- imputation_file %>%
  select(
    ebird_COMMON_NAME, max_distance, max_brightness,
    mass_log10, flock_size_log10, IUCN_ordinal
  )

joined_dat <- dat %>%
  left_join(dat_imp, by = "ebird_COMMON_NAME") %>%
  left_join(detectability, by = "ebird_COMMON_NAME") %>%
  mutate(
    density = 10^(intercept + slope * log10_mean_abund),
    log10_density = log10(density)
  )

####################################
# 7. UNCERTAINTY FUNCTIONS
####################################

predict_se <- function(beta_sim, x) {
  xb <- cbind(1, x) %*% t(beta_sim)
  qs <- apply(xb, 1, quantile, probs = c(0.5, 0.8413448))
  qs[2, ] - qs[1, ]
}

brms_mod <- readRDS("Intermediate data files/brms_mod.rds")

coef_names <- colnames(vcov(brms_mod))
intercept_name <- coef_names[coef_names == "Intercept"]
slope_name <- coef_names[coef_names != "Intercept"]
stopifnot(length(slope_name) == 1)

post <- as_draws_df(brms_mod)

simulate_gaussian <- function(n = 10000) {
  vcv <- vcov(brms_mod)[
    c(intercept_name, slope_name),
    c(intercept_name, slope_name)
  ]
  mu <- fixef(brms_mod)[
    c(intercept_name, slope_name),
    "Estimate"
  ]
  rmvnorm(n, mu, vcv)
}

####################################
# 8. RUN UNCERTAINTY PIPELINE
####################################

training_dat <- joined_dat %>% filter(!is.na(log10_density))
data_for_imputation <- joined_dat %>% filter(is.na(log10_density))

grid_area <- readRDS("Data/Spatial_data/grid_percent_land.RDS") %>%
  select(ID, percentage) %>%
  rename(grid_id = ID) %>%
  st_set_geometry(NULL)

UNCERTAINTY_METHODS <- list(
  gaussian = function(d) {
    simulate_gaussian()
  }
)


for (method in names(UNCERTAINTY_METHODS)) {
  
  message("Running uncertainty method: ", method)
  
  SEs <- bind_rows(lapply(unique(training_dat$ebird_COMMON_NAME), function(sp) {
    
    spp <- training_dat %>% filter(ebird_COMMON_NAME == sp)
    
    bind_rows(lapply(unique(spp$grid_id), function(g) {
      
      d <- spp %>% filter(grid_id == g)
      beta_sim <- UNCERTAINTY_METHODS[[method]](d)
      
      data.frame(
        ebird_COMMON_NAME = sp,
        grid_id = g,
        se = predict_se(beta_sim, d$log10_mean_abund)
      )
    }))
  }))
  
  out <- bind_rows(
    training_dat %>% left_join(SEs, by = c("ebird_COMMON_NAME", "grid_id")),
    data_for_imputation %>% mutate(se = NA)
  ) %>%
    left_join(grid_area, by = "grid_id")
  
  attr(out, "uncertainty_method") <- method
  attr(out, "generated_on") <- Sys.time()
  
  saveRDS(
    out,
    paste0("Intermediate data files/file_for_imputation_", method, ".RDS")
  )
}