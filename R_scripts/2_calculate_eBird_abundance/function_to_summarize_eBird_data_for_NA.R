### This is an R script for a minimum working example
### in order to get relative abundance for a species
### across many states and the necessary bcrs
### I have already chunked all the effort data for each state from eBird
### and species-specific data so should have everything we need to build models
### will source this function and then loop through the species
### in order to calculate a RDS for each species
library(readr)
## read in cleaned pop estimates
pop_estimates <- read_csv("Data/PIF estimates/cleaned_pop_estimates.csv")

# large overall function for a given species
summarize_abundances_for_a_species <- function(species) {

  message(paste0("Analyzing ", species))
  
  ## packages
  require(readr)
  require(dplyr)
  require(tidyr)
  require(mgcv)
  require(lubridate)

  # subset and manipulate pop_estimates for those
  # locations and the species that we have data for
  pop_estimates2 <- pop_estimates %>%
    dplyr::filter(English_Name == species) %>%
    mutate(COUNTRY = gsub("USA", "US", Country)) %>%
    mutate(COUNTRY = gsub("CAN", "CA", COUNTRY)) %>%
    unite(STATE_CODE_BCR, COUNTRY, PROV_STATE, BCR, sep="-", remove=FALSE)

  # read in bird data
  bird_dat <- readRDS(paste0("eBird data/species_dat/", gsub(" ", "_", species), ".RDS"))


  # get state codes from where we know population estimates
  # these state codes are where we'll then generate eBird estimates from
  state_code_bcr_list <- pop_estimates2 %>%
    mutate(COUNTRY = gsub("USA", "US", Country)) %>%
    mutate(COUNTRY = gsub("CAN", "CA", COUNTRY)) %>%
    unite(STATE_CODE_BCR, COUNTRY, PROV_STATE, BCR, sep="-", remove=FALSE) %>%
    unite(STATE_CODE, COUNTRY, PROV_STATE, sep="-") %>%
    dplyr::select(STATE_CODE_BCR, STATE_CODE) %>%
    distinct(STATE_CODE_BCR, .keep_all = TRUE)

  # filter bird data to include only data from these state codes
  # we aren't interested in generating estimates from other
  # localities yet! Only those we have training data for
  bird_dat2 <- bird_dat %>%
    unite(STATE_CODE_BCR, STATE_CODE, BCR_CODE, sep="-", remove=FALSE) %>%
    dplyr::filter(STATE_CODE_BCR %in% state_code_bcr_list$STATE_CODE_BCR)

  # now for a given state code
  # read in the effort data and run a model
  # this is a nested function that applies the below functions for each state_code
  # that the species exists within
run_function_for_each_state <- function(state) {
  
  state_effort <- readRDS(paste0("eBird data/effort_dat/", state, ".RDS")) %>%
    unite(STATE_CODE_BCR, STATE_CODE, BCR_CODE, sep="-", remove=FALSE) 
  
  # this is a mid-level function
  # that splits the state up by bcr and runs the model
  run_function_by_state_code_bcr <- function(state_code_bcr){
    
    message(paste0("Analyzing ", state_code_bcr))
    
    analysis_dat <- suppressMessages(state_effort %>%
      dplyr::filter(STATE_CODE_BCR == state_code_bcr) %>%
      left_join(., bird_dat2) %>%
      replace_na(list(OBSERVATION_COUNT = 0)) %>%
      dplyr::filter(OBSERVATION_COUNT != "X") %>%
      mutate(OBSERVATION_COUNT = as.numeric(as.character(OBSERVATION_COUNT))) %>%
      mutate(month=month(OBSERVATION_DATE, label=TRUE)) %>%
      mutate(month=as.character(as.factor(month)))) %>%
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
        mutate(STATE_CODE_BCR = state_code_bcr) %>%
        mutate(COMMON_NAME = species) %>%
        mutate(MONTH = month_of_year) %>%
        mutate(STATE_CODE = state) %>%
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
    
    bottom_level_df <- bind_rows(results_list_by_month)
    
    return(bottom_level_df)
    
  }
  
  # run the mid level function by each bcr code state combination
  state_code_bcr_list2 <- state_code_bcr_list %>%
    dplyr::filter(STATE_CODE == state) %>%
    .$STATE_CODE_BCR
  
  results_list_by_state_code_bcr <- lapply(state_code_bcr_list2, function(x){run_function_by_state_code_bcr(x)})
  
  state_code_bcr_df <- bind_rows(results_list_by_state_code_bcr)
  
  return(state_code_bcr_df)
}


list_of_states <- unique(state_code_bcr_list$STATE_CODE)

results_list_by_state <- lapply(list_of_states, function(x){run_function_for_each_state(x)})


top_level_df <- bind_rows(results_list_by_state)


## Create a 'final dataset' for a given species
## combining with the population estimates from PIF
## and save this as an RDS
## by saving each species, it means it is less likely to have issues and can then
## not have to loop through all species every time

final_dataset <- pop_estimates2 %>%
  dplyr::select(-COUNTRY, -Sequence_AOS_59) %>%
  rename(COMMON_NAME = English_Name) %>%
  rename(BCR_CODE = BCR) %>%
  left_join(., top_level_df, by=c("COMMON_NAME", "STATE_CODE_BCR"))
  

saveRDS(final_dataset, file=paste0("Data/analyzed_eBird_abundances/NA_", gsub(" ", "_", species), ".RDS"))

}




