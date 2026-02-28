### This is an R script for a minimum working example
### in order to get relative abundance for a species
### across many states and the necessary bcrs
### I have already chunked all the effort data for each state from eBird
### and species-specific data so should have everything we need to build models
### will source this function and then loop through the species
### in order to calculate a RDS for each species
library(readr)
## read in cleaned pop estimates
pop_estimates <- read_csv("Data/Species range maps/ASSESSED_potential_birdlife_species_for_inclusion.csv")

# large overall function for a given species
summarize_abundances_for_a_species <- function(species) {
  
  message(paste0("Analyzing ", species))
  
  ## packages
  require(readr)
  require(dplyr)
  require(tidyr)
  require(lubridate)
  
  # read in bird data
  bird_dat <- readRDS(paste0("eBird data/species_dat/", gsub(" ", "_", species), "_Birdlife_dat.RDS"))

  # get a list of sampling event identifiers that the species occurs on
  # this needs to be done separately because Xs can be submitted for the
  # species of interest
  # and I'll get rid of those
  lists_to_get_rid_of <- bird_dat %>%
    dplyr::filter(CATEGORY %in% c("species", "issf", "domestic")) %>%
    dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
    dplyr::filter(COMMON_NAME == species) %>%
    .$SAMPLING_EVENT_IDENTIFIER
  
  # get all the lists with the species dat
  # but get rid of any lists which had an X for the species of interest
  lists_with_species_dat <- bird_dat %>%
    dplyr::filter(CATEGORY %in% c("species", "issf", "domestic")) %>%
    dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
    dplyr::filter(COMMON_NAME == species) %>%
    mutate(OBSERVATION_COUNT = as.numeric(as.character(OBSERVATION_COUNT))) %>%
    dplyr::filter(complete.cases(OBSERVATION_COUNT)) %>%
    dplyr::select(-COMMON_NAME, -SCIENTIFIC_NAME)
  
  # get the rest of the lists (i.e., effort dat)
  all_other_lists <- bird_dat %>%
    dplyr::filter(CATEGORY %in% c("species", "issf", "domestic")) %>%
    dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
    dplyr::filter(!SAMPLING_EVENT_IDENTIFIER %in% lists_to_get_rid_of) %>%
    dplyr::select(SAMPLING_EVENT_IDENTIFIER, STATE_CODE, LOCALITY_ID, LATITUDE, LONGITUDE,
                  OBSERVATION_DATE, TIME_OBSERVATIONS_STARTED, OBSERVER_ID,
                  PROTOCOL_TYPE, DURATION_MINUTES, EFFORT_DISTANCE_KM, EFFORT_AREA_HA,
                  NUMBER_OBSERVERS, CATEGORY, GROUP_IDENTIFIER) %>%
    mutate(OBSERVATION_COUNT = 0) %>%
    distinct()
  
  # combine the two dataframes into an 'analysis_dat' dataframe
  analysis_dat <- bind_rows(lists_with_species_dat,
                            all_other_lists) %>%
    mutate(month=month(OBSERVATION_DATE, label=TRUE)) %>%
    mutate(month=as.character(as.factor(month))) %>%
    replace_na(list(EFFORT_DISTANCE_KM=0))
  
  # this is the bottom level function that runs the model by each month
  # and summarizes the data observations for that month
  run_model_by_month_function <- function(month_of_year) {
    
    message(paste0("Analyzing ", month_of_year))
    
    dat_for_model <- analysis_dat %>%
      dplyr::filter(month==month_of_year) %>%
      tibble::rownames_to_column(var="rowID") %>%
      mutate(GROUP_IDENTIFIER.2 = case_when(
        is.na(GROUP_IDENTIFIER) == TRUE ~ rowID,
        is.na(GROUP_IDENTIFIER) == FALSE ~ GROUP_IDENTIFIER
      )) %>%
      distinct(GROUP_IDENTIFIER.2, .keep_all=TRUE)
    
    glm_mod <- glm(OBSERVATION_COUNT ~ DURATION_MINUTES + EFFORT_DISTANCE_KM, family="poisson",
                   data=dat_for_model)
    
    pd <- data.frame(DURATION_MINUTES = 60, EFFORT_DISTANCE_KM = 1)
    
    answer <- as.data.frame(broom::augment(x=glm_mod, newdata = pd, type.predict = "response", se_fit=TRUE)) %>%
      rename(model_fit=.fitted) %>%
      rename(model_se=.se.fit) %>%
      mutate(COMMON_NAME = species) %>%
      mutate(MONTH = month_of_year) %>%
      mutate(sum_of_abundance=sum(dat_for_model$OBSERVATION_COUNT)) %>%
      mutate(sd_of_abundance=sd(dat_for_model$OBSERVATION_COUNT)) %>%
      mutate(mean_of_abundance=mean(dat_for_model$OBSERVATION_COUNT)) %>%
      mutate(median_of_abundance=median(dat_for_model$OBSERVATION_COUNT)) %>%
      mutate(total_effort_time=sum(dat_for_model$DURATION_MINUTES, na.rm=TRUE)) %>%
      mutate(total_effort_distance=sum(dat_for_model$EFFORT_DISTANCE_KM, na.rm=TRUE)) %>%
      mutate(mean_effort_time=mean(dat_for_model$DURATION_MINUTES, na.rm=TRUE)) %>%
      mutate(mean_effort_distance=mean(dat_for_model$EFFORT_DISTANCE_KM, na.rm=TRUE)) %>%
      mutate(number_obs_in_model = nrow(dat_for_model)) %>%
      mutate(model_converged = glm_mod$converged)
    
  }
  
  # run the above portion of the model for every month possible
  # in the data
  months_available <- unique(analysis_dat$month)
  
  
  results_list_by_month <- lapply(months_available, function(x){run_model_by_month_function(x)})
  
  results_df <- bind_rows(results_list_by_month)
  
  
  ## Create a 'final dataset' for a given species
  ## combining with the population estimates from the british birds estimate
  ## and save this as an RDS
  ## by saving each species, it means it is less likely to have issues and can then
  ## not have to loop through all species every time
  
  final_dataset <- pop_estimates %>%
    right_join(., results_df, by="COMMON_NAME")
  
  print(bird_dat)
  saveRDS(final_dataset, file=paste0("Data/analyzed_eBird_abundances/Birdlife_", gsub(" ", "_", species), ".RDS"))
  
}
