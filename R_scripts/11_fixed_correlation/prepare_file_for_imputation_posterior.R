###############################################################
# 0. Packages
###############################################################

library(brms)
library(dplyr)
library(purrr)

###############################################################
# 1. Load brms model
###############################################################

fit <- readRDS("Intermediate data files/brms_mod.rds")

re <- ranef(fit, summary = FALSE)

species_array <- re$species  # draws × species × 2
species_names <- dimnames(species_array)[[2]]

intercept_draws <- species_array[, , "Intercept"]
slope_draws     <- species_array[, , "meabundance_log10abundance_SD_log10"]

clements_clean <- read_csv("Data/clements_clean.csv")

body_size <- read_csv("Data/body size data/cleaned_body_size_data.csv") %>%
  select(3, 7:9)

color <- read_csv("Data/color data/color_data_by_species.csv")

IUCN <- read_csv("Data/IUCN categories/cleaned_IUCN_data.csv")

flock_size <- readRDS("Data/Flock size/flock_size_per_month_for_each_species.RDS") %>%
  group_by(COMMON_NAME) %>%
  summarise(mean_max_flock_size = mean(max_abund)) %>%
  rename(ebird_COMMON_NAME = COMMON_NAME)

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
    flock_size_log10 = log10(mean_max_flock_size)
  ) %>%
  ungroup()

months <- c("january","february","march","april","may","june",
            "july","august","september","october","november","december")

combined_dat <- map_df(months, ~ readRDS(
  paste0("Intermediate data files/abundance_grid_results/", .x, ".RDS")
))

temporal_summary <- combined_dat %>%
  filter(mean_abund > 0) %>%
  group_by(grid_id, ebird_COMMON_NAME) %>%
  summarise(
    mean_abund = mean(mean_abund),
    number_months = n_distinct(MONTH),
    mean_species_checklists = mean(obs_with_abund),
    total_species_checklists = sum(obs_with_abund),
    .groups = "drop"
  ) %>%
  left_join(
    combined_dat %>%
      filter(mean_abund > 0) %>%
      select(grid_id, total_checklists) %>%
      distinct() %>%
      group_by(grid_id) %>%
      summarise(
        mean_grid_checklists = mean(total_checklists),
        total_grid_checklists = sum(total_checklists),
        .groups = "drop"
      ),
    by = "grid_id"
  ) %>%
  left_join(
    combined_dat %>%
      select(grid_id, area_square_miles) %>%
      distinct(),
    by = "grid_id"
  )

rm(combined_dat)

joined_dat <- temporal_summary %>%
  left_join(imputation_file, by = "ebird_COMMON_NAME") %>%
  mutate(
    species = ebird_COMMON_NAME,
    abundance_log10 = log10(mean_abund),
    log10_mean_abund = abundance_log10
  )

getSE_fast <- function(sp){
  
  spp_data <- joined_dat[joined_dat$species == sp, ]
  
  if(nrow(spp_data) == 0) return(NULL)
  
  sp_index <- match(sp, species_names)
  
  int_draws <- intercept_draws[, sp_index]
  slo_draws <- slope_draws[, sp_index]
  
  x <- spp_data$abundance_log10
  
  eta <- outer(int_draws, rep(1, length(x))) +
    outer(slo_draws, x)
  
  data.frame(
    grid_id = spp_data$grid_id,
    species = sp,
    log10_density = apply(eta, 2, median),
    se = apply(eta, 2, sd)
  )
}

species_list <- intersect(unique(joined_dat$species),
                          species_names)

SEs_1 <- do.call(
  rbind,
  lapply(species_list, getSE_fast)
)

joined_dat_final <- joined_dat %>%
  left_join(SEs_1,
            by = c("grid_id", "species"))

saveRDS(joined_dat_final,
        "Intermediate data files/file_for_imputation_posterior.RDS")
