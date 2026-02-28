######## This is an R script used to
######## do taxonomic matching between four different sources of data
######## that will be used in our analyses throughout
######## 1) eBird/clements taxonomy
######## 2) Jetz's phylogenetic tree
######## 3) Body size data
######## 4) Bird color data
######## 5) IUCN status


# packages
library(readr)
library(dplyr)
library(tidyr)

clements_list <- read_csv("Data/Clements-Checklist-v2018-August-2018.csv") %>%
  dplyr::filter(category == "species") %>%
  rename(COMMON_NAME = `English name`) %>%
  rename(SCIENTIFIC_NAME = `scientific name`) %>%
  dplyr::select(COMMON_NAME, SCIENTIFIC_NAME, order, family)


phylo_tree <- read_csv("Data/phylo/2012-03-04206D-master_taxonomy.csv") %>%
  rename(SCIENTIFIC_NAME = Scientific) %>%
  rename(COMMON_NAME = English)


body_size <- read_csv("Data/body size data/Amniote_Database_Aug_2015.csv") %>%
  rename(COMMON_NAME = common_name) %>%
  unite(SCIENTIFIC_NAME, genus, species, sep=" ")


color <- read_csv("Data/color data/Dale_2015_data_allpoints.csv")

IUCN <- read_csv("Data/IUCN categories/HBW-BirdLife_Checklist_Version_3.csv") %>%
  dplyr::select(4,5,8) %>%
  rename(COMMON_NAME=`Common name`) %>%
  rename(SCIENTIFIC_NAME=`Scientific name`) %>%
  rename(IUCN_category=`2018 IUCN Red List category`) %>%
  replace_na(list(COMMON_NAME = "SUBSPECIES")) %>%
  dplyr::filter(COMMON_NAME != "SUBSPECIES")


### Everything needs to match back to the phylogenetic tree
### AS in, this will be our destination
### so, I'll go one at a time through the tree and see how many don't match up
### with each of the datasets

# test color first
color_join <- color %>%
  dplyr::select(TipLabel) %>%
  distinct() %>%
  mutate(color_data = "Yes") %>%
  left_join(., phylo_tree, by="TipLabel")

# looks like color ALL match, based on TipLabel 

# now test body size
body_size_join_scientific <- body_size %>%
  dplyr::select(SCIENTIFIC_NAME) %>%
  distinct() %>%
  mutate(body_size_data = "Yes") %>%
  left_join(., phylo_tree, by="SCIENTIFIC_NAME") %>%
  replace_na(list(TipLabel="MISSING"))

missing_body_size_scientific <- body_size_join_scientific %>%
  dplyr::filter(TipLabel == "MISSING") %>%
  dplyr::select(SCIENTIFIC_NAME) %>%
  dplyr::left_join(., body_size, by="SCIENTIFIC_NAME")

body_size_scientific_matched <- body_size_join_scientific %>%
  dplyr::filter(TipLabel != "MISSING") %>%
  dplyr::select(-body_size_data) %>%
  dplyr::select(SCIENTIFIC_NAME, TipLabel)


# test the missing portion based on common name
body_size_join_common <- missing_body_size_scientific %>%
  dplyr::select(COMMON_NAME) %>%
  distinct() %>%
  mutate(body_size_data = "Yes") %>%
  left_join(., phylo_tree, by="COMMON_NAME") %>%
  replace_na(list(TipLabel="MISSING"))

missing_body_size_common <- body_size_join_common %>%
  dplyr::filter(TipLabel == "MISSING") %>%
  dplyr::select(COMMON_NAME) %>%
  dplyr::left_join(., body_size, by="COMMON_NAME")

body_size_common_matched <- body_size_join_common %>%
  dplyr::filter(TipLabel != "MISSING") %>%
  dplyr::select(-body_size_data) %>%
  dplyr::select(COMMON_NAME, TipLabel)

body_size_cleaned <- body_size %>%
  dplyr::select(1:4, 6, 10, 22, 23) %>%
  left_join(., body_size_scientific_matched, by="SCIENTIFIC_NAME") %>%
  left_join(., body_size_common_matched, by="COMMON_NAME") %>%
  mutate(TipLabel = ifelse(is.na(TipLabel.x) == "FALSE", TipLabel.x, TipLabel.y)) %>%
  dplyr::select(-TipLabel.x, -TipLabel.y) %>%
  filter(complete.cases(TipLabel)) %>%
  dplyr::select(COMMON_NAME, SCIENTIFIC_NAME, TipLabel, 1:3, 6:8)


write_csv(body_size_cleaned, "Data/body size data/cleaned_body_size_data.csv")


# now repeat the above, but for IUCN categories
IUCN_join_scientific <- IUCN %>%
  dplyr::select(SCIENTIFIC_NAME) %>%
  distinct() %>%
  mutate(IUCN_data = "Yes") %>%
  left_join(., phylo_tree, by="SCIENTIFIC_NAME") %>%
  replace_na(list(TipLabel="MISSING"))

missing_IUCN_scientific <- IUCN_join_scientific %>%
  dplyr::filter(TipLabel == "MISSING") %>%
  dplyr::select(SCIENTIFIC_NAME) %>%
  dplyr::left_join(., IUCN, by="SCIENTIFIC_NAME")

IUCN_scientific_matched <- IUCN_join_scientific %>%
  dplyr::filter(TipLabel != "MISSING") %>%
  dplyr::select(-IUCN_data) %>%
  dplyr::select(SCIENTIFIC_NAME, TipLabel)


# test the missing portion based on common name
IUCN_join_common <- missing_IUCN_scientific %>%
  dplyr::select(COMMON_NAME) %>%
  distinct() %>%
  mutate(IUCN_data = "Yes") %>%
  left_join(., phylo_tree, by="COMMON_NAME") %>%
  replace_na(list(TipLabel="MISSING"))

missing_IUCN_common <- IUCN_join_common %>%
  dplyr::filter(TipLabel == "MISSING") %>%
  dplyr::select(COMMON_NAME) %>%
  dplyr::left_join(., IUCN, by="COMMON_NAME")

IUCN_common_matched <- IUCN_join_common %>%
  dplyr::filter(TipLabel != "MISSING") %>%
  dplyr::select(-IUCN_data) %>%
  dplyr::select(COMMON_NAME, TipLabel)

IUCN_cleaned <- IUCN %>%
  left_join(., IUCN_scientific_matched, by="SCIENTIFIC_NAME") %>%
  left_join(., IUCN_common_matched, by="COMMON_NAME") %>%
  mutate(TipLabel = ifelse(is.na(TipLabel.x) == "FALSE", TipLabel.x, TipLabel.y)) %>%
  dplyr::select(-TipLabel.x, -TipLabel.y) %>%
  filter(complete.cases(TipLabel))


write_csv(IUCN_cleaned, "Data/IUCN categories/cleaned_IUCN_data.csv")





### Now for the tricky one! eBird data
### but do it from tree to eBird
### and then see how many are missing
### and then fix those...
### Really, I think it will be okay to
### get a one-way join.
### As long as all 9,993 tips of the tree are associated with
### one eBird/clements species, then we can do something simple
### where we assign 'missing' clements species, the slope based on
### the other slopes in the genus
ebird_join_scientific <- clements_list %>%
  dplyr::select(SCIENTIFIC_NAME) %>%
  distinct() %>%
  mutate(ebird = "Yes") %>%
  left_join(., phylo_tree, by="SCIENTIFIC_NAME") %>%
  replace_na(list(TipLabel="MISSING"))

missing_ebird_scientific <- ebird_join_scientific %>%
  dplyr::filter(TipLabel == "MISSING") %>%
  dplyr::select(SCIENTIFIC_NAME) %>%
  dplyr::left_join(., clements_list, by="SCIENTIFIC_NAME")

ebird_scientific_matched <- ebird_join_scientific %>%
  dplyr::filter(TipLabel != "MISSING") %>%
  dplyr::select(-ebird) %>%
  dplyr::select(SCIENTIFIC_NAME, TipLabel) %>%
  rename(ebird_SCIENTIFIC_NAME=SCIENTIFIC_NAME)


# test the missing portion based on common name
ebird_join_common <- missing_ebird_scientific %>%
  dplyr::select(COMMON_NAME) %>%
  distinct() %>%
  mutate(ebird = "Yes") %>%
  left_join(., phylo_tree, by="COMMON_NAME") %>%
  replace_na(list(TipLabel="MISSING"))

missing_ebird_common <- ebird_join_common %>%
  dplyr::filter(TipLabel == "MISSING") %>%
  dplyr::select(COMMON_NAME) %>%
  dplyr::left_join(., clements_list, by="COMMON_NAME")

ebird_common_matched <- ebird_join_common %>%
  dplyr::filter(TipLabel != "MISSING") %>%
  dplyr::select(-ebird) %>%
  dplyr::select(COMMON_NAME, TipLabel) %>%
  rename(ebird_COMMON_NAME = COMMON_NAME)


matched_birds <- ebird_common_matched %>%
  dplyr::select(TipLabel) %>%
  mutate(matched="YES") %>%
  bind_rows(., dplyr::select(mutate(ebird_scientific_matched, matched="YES"), TipLabel, matched)) %>%
  distinct()


non_matched_tips <- phylo_tree %>%
  left_join(., matched_birds, by="TipLabel") %>%
  replace_na(list(matched="No")) %>%
  dplyr::filter(matched == "No") %>%
  rename(SCIENTIFIC_NAME_tree=SCIENTIFIC_NAME) %>%
  rename(COMMON_NAME_tree=COMMON_NAME)

## Export csv for manual matching outside of R
write_csv(non_matched_tips, "Data/phylo/non_matched_tips_for_manual_matching_with_clements.csv")

## read in file with manually matched different birds
manually_matched_tips <- read_csv("Data/phylo/FIXED_manual_matched_tips_with_clements.csv") %>%
  dplyr::select(ebird_SCIENTIFIC_NAME, TipLabel)
  
  
clements_clean <- clements_list %>%
  rename(ebird_SCIENTIFIC_NAME = SCIENTIFIC_NAME) %>%
  rename(ebird_COMMON_NAME = COMMON_NAME) %>%
  left_join(., ebird_scientific_matched, by="ebird_SCIENTIFIC_NAME") %>%
  left_join(., ebird_common_matched, by="ebird_COMMON_NAME") %>%
  mutate(TipLabel = ifelse(is.na(TipLabel.x) == "FALSE", TipLabel.x, TipLabel.y)) %>%
  dplyr::select(-TipLabel.x, -TipLabel.y) %>%
  left_join(., manually_matched_tips, by="ebird_SCIENTIFIC_NAME") %>%
  mutate(TipLabel = ifelse(is.na(TipLabel.y) == "TRUE", TipLabel.x, TipLabel.y)) %>%
  dplyr::select(-TipLabel.x, -TipLabel.y)
  
write_csv(clements_clean, "Data/clements_clean.csv")




