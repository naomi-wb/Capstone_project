# this is an R script to get the mean abundance of all listed species in each
# grid, where the mean abundance is the mean of abundance across all checklists

## packages
library(readr)
library(DBI)
library(bigrquery)
library(dbplyr)
library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(scales)
library(patchwork)


# create connection with online database
con <- DBI::dbConnect(bigrquery::bigquery(),
                      dataset= "ebird",
                      project="ebird-database",
                      billing="ebird-database")

# create ebird table
ebird <- tbl(con, 'ebird_qa')

# create locality_id and grid_id table
grids <- tbl(con, 'grid_5_joined_with_sites')

# first query database for the grid/species/month combination for the total abundance
# on all lists which did not have an X associated with them
# this removes the 'x' lists from contention
grids1 <- ebird %>%
  dplyr::select(COMMON_NAME, SCIENTIFIC_NAME, SAMPLING_EVENT_IDENTIFIER, LOCALITY_ID,
                OBSERVATION_COUNT, OBSERVATION_DATE, SCIENTIFIC_NAME) %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  left_join(., grids, by="LOCALITY_ID") %>%
  mutate(MONTH=month(OBSERVATION_DATE)) %>%
  dplyr::filter(OBSERVATION_COUNT != "X") %>%
  mutate(OBSERVATION_COUNT=as.numeric(as.character(OBSERVATION_COUNT))) %>%
  group_by(MONTH, grid_id, COMMON_NAME, SCIENTIFIC_NAME) %>%
  summarize(total_sum=sum(OBSERVATION_COUNT), 
            obs_with_abund=n()) %>%
  collect(n=Inf)

# then query the database and for each grid/species/month
# get the number of X's in that category
grids2 <- ebird %>%
  dplyr::select(COMMON_NAME, SCIENTIFIC_NAME, SAMPLING_EVENT_IDENTIFIER, LOCALITY_ID,
                OBSERVATION_COUNT, OBSERVATION_DATE, SCIENTIFIC_NAME) %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  left_join(., grids, by="LOCALITY_ID") %>%
  mutate(MONTH=month(OBSERVATION_DATE)) %>%
  dplyr::filter(OBSERVATION_COUNT == "X") %>%
  group_by(MONTH, grid_id, COMMON_NAME, SCIENTIFIC_NAME) %>%
  summarize(obs_with_x=n()) %>%
  collect(n=Inf)

# now get the total number of lists
# in each grid/month combination
# which will provide us with how many '0s' there were when joined with
# the above two datasets
grid_effort <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, LOCALITY_ID, OBSERVATION_DATE) %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>% 
  left_join(., grids, by="LOCALITY_ID") %>%
  mutate(MONTH=month(OBSERVATION_DATE)) %>%
  dplyr::select(grid_id, MONTH, SAMPLING_EVENT_IDENTIFIER) %>%
  distinct() %>%
  group_by(grid_id, MONTH) %>%
  summarize(total_checklists=n()) %>%
  collect(n=Inf)

length(unique(grid_effort$grid_id))

# plot a histogram of this effort per grid cell, per month
grid_effort_plot <- ggplot(grid_effort, aes(x=total_checklists))+
  geom_histogram(fill="lightblue", color="black", bins=50)+
  scale_x_log10(labels=comma)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Number of checklists in a grid")+
  ylab("Number of grids")+
  ggtitle("A")

grid_effort_plot

ggsave("Figures/Supplementary_figures/histogram_of_checklists_per_grid_month.png", 
       width=7.8, height=5.9, units="in")

# now remove grids with <50 eBird checklists
# and replot the histogram
grid_effort_trimmed <- grid_effort %>%
  dplyr::filter(total_checklists>=50)

grid_effort_trimmed_plot <- ggplot(grid_effort_trimmed, aes(x=total_checklists))+
  geom_histogram(fill="lightblue", color="black", bins=50)+
  scale_x_log10(labels=comma)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Number of checklists in a grid")+
  ylab("Number of grids")+
  ggtitle("B")

grid_effort_trimmed_plot

ggsave("Figures/Supplementary_figures/histogram_of_checklists_per_grid_month_trimmed_at_50.png", 
       width=7.8, height=5.9, units="in")

length(unique(grid_effort_trimmed$grid_id))

grid_effort_plot + grid_effort_trimmed_plot + plot_layout(ncol=1)

ggsave("Figures/Supplementary_figures/histograms_of_checklists_per_grid.png", 
       width=5.8, height=7.9, units="in")

# now need to work with these three dataframes
# and this will help to provide the final 'answers'
# that we need for summarizing abundance in the world
grids3 <- grids1 %>% 
  left_join(., grids2) %>%
  replace_na(list(obs_with_x=0)) %>%
  left_join(., grid_effort) %>%
  mutate(mean_abund=(total_sum/(total_checklists-obs_with_x))) %>%
  dplyr::filter(mean_abund != "Inf") # this happens when total checklists and obs with x are equal, so remove

# now read in the original grids geojson
# to see about calculating the 'area' of each grid
grids_5_degree <- st_read("Data/Spatial_data/grid_5_degree.geojson")

length(unique(grids_5_degree$ID))

# get a dataframe of grid_id and the area of the grid
grid_area <- data.frame(grid_id=grids_5_degree$ID,
                        area=st_area(grids_5_degree)) %>%
  mutate(area_square_miles=as.numeric(area/2.59e+6)) %>%
  dplyr::select(grid_id, area_square_miles)

# now we need to remove all the weird 'categories'
# in eBird
# and also zero fill every grid with the species that are found elsewhere
# will use our clean clements dataset for this
clements <- read_csv("Data/clements_clean.csv") %>%
  group_by(ebird_COMMON_NAME, ebird_SCIENTIFIC_NAME) %>%
  slice(1) %>%
  ungroup()

# first remove any "COMMON_NAME" from grids3
# that doesn't match with a CLEMENTS value from that csv
# original unique common names
length(unique(grids3$COMMON_NAME))

grids4 <- grids3 %>%
  dplyr::filter(COMMON_NAME %in% clements$ebird_COMMON_NAME)

# common names after filtering
length(unique(grids4$COMMON_NAME)) # this makes a lot more sense!

# now join with grid_area
grids5 <- grids4 %>%
  left_join(., grid_area)

# now make a plot of all grids with checklists
# and then a plot of grids after trimming by checklists
# download some natural earth data
countries <- ne_countries(returnclass = "sf") %>%
  st_geometry()

graticules <- ne_download(type = "graticules_15", category = "physical",
                          returnclass = "sf") %>%
  st_geometry()

bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>%
  st_geometry()


# make a plot of the grids over the world
# for both possible scenarios to visualize it
# which may come in handy for supplementary material
all_possible_grids <- grids_5_degree %>%
  rename(grid_id=ID) %>%
  left_join(., grid_effort %>%
              ungroup() %>%
              dplyr::select(grid_id) %>%
              distinct() %>%
              mutate(Sampled="True")) %>%
  replace_na(list(Sampled="False"))

trimmed_grids <- grids_5_degree %>%
  rename(grid_id=ID) %>%
  left_join(., grid_effort %>%
              dplyr::filter(total_checklists>=50) %>%
              ungroup() %>%
              dplyr::select(grid_id) %>%
              distinct() %>%
              mutate(Sampled="True")) %>%
  replace_na(list(Sampled="False"))

# write out the grids which were sampled
# as a geojson for use later
st_write(trimmed_grids, "Data/Spatial_data/grids_sampled_by_50_checklists.geojson")

all_grids_plot <- ggplot()+
  geom_sf(data = bb, col = "grey20", fill = "transparent")+
  geom_sf(data = countries, fill = "black", col = "grey40", lwd = 0.3)+
  geom_sf(data = all_possible_grids, aes(fill=Sampled), alpha=0.8)+
  theme_minimal()+
  theme(axis.text = element_blank())+
  scale_fill_brewer(palette="Set2")+
  theme(legend.position="bottom")+
  ggtitle("A   All grids with at least one eBird checklist in at least one month")

all_grids_plot

ggsave("Figures/Supplementary_figures/all_grids_with_at_least_one_checklist.png", 
       width=8.3, height=5.9, units="in")

trimmed_grids_plot <- ggplot()+
  geom_sf(data = bb, col = "grey20", fill = "transparent")+
  geom_sf(data = countries, fill = "black", col = "grey40", lwd = 0.3)+
  geom_sf(data = trimmed_grids, aes(fill=Sampled), alpha=0.8)+
  theme_minimal()+
  theme(axis.text = element_blank())+
  scale_fill_brewer(palette="Set2")+
  theme(legend.position="bottom")+
  ggtitle("B   All grids with at least 50 eBird checklists in at least one month")

trimmed_grids_plot

ggsave("Figures/Supplementary_figures/all_grids_with_at_least_fifty_checklists.png", 
       width=8.3, height=5.9, units="in")

all_grids_plot + trimmed_grids_plot + plot_layout(ncol=1)

ggsave("Figures/Supplementary_figures/spatial_plot_of_grids.png", 
       width=5.8, height=7.9, units="in")

# now trim the entire grids dataset
grids6 <- grids5 %>%
  dplyr::filter(total_checklists>=50)


# now split by month and save out the dataframes
# as individual RDSs
Jan <- grids6 %>%
  dplyr::filter(MONTH==1)
saveRDS(Jan, "Intermediate data files/global_mean_abundance_from_ebird/january.RDS")

Feb <- grids6 %>%
  dplyr::filter(MONTH==2)
saveRDS(Feb, "Intermediate data files/global_mean_abundance_from_ebird/february.RDS")

Mar <- grids6 %>%
  dplyr::filter(MONTH==3)
saveRDS(Mar, "Intermediate data files/global_mean_abundance_from_ebird/march.RDS")

Apr <- grids6 %>%
  dplyr::filter(MONTH==4)
saveRDS(Apr, "Intermediate data files/global_mean_abundance_from_ebird/april.RDS")

May <- grids6 %>%
  dplyr::filter(MONTH==5)
saveRDS(May, "Intermediate data files/global_mean_abundance_from_ebird/may.RDS")

Jun <- grids6 %>%
  dplyr::filter(MONTH==6)
saveRDS(Jun, "Intermediate data files/global_mean_abundance_from_ebird/june.RDS")

Jul <- grids6 %>%
  dplyr::filter(MONTH==7)
saveRDS(Jul, "Intermediate data files/global_mean_abundance_from_ebird/july.RDS")

Aug <- grids6 %>%
  dplyr::filter(MONTH==8)
saveRDS(Aug, "Intermediate data files/global_mean_abundance_from_ebird/august.RDS")

Sep <- grids6 %>%
  dplyr::filter(MONTH==9)
saveRDS(Sep, "Intermediate data files/global_mean_abundance_from_ebird/september.RDS")

Oct <- grids6 %>%
  dplyr::filter(MONTH==10)
saveRDS(Oct, "Intermediate data files/global_mean_abundance_from_ebird/october.RDS")

Nov <- grids6 %>%
  dplyr::filter(MONTH==11)
saveRDS(Nov, "Intermediate data files/global_mean_abundance_from_ebird/november.RDS")

Dec <- grids6 %>%
  dplyr::filter(MONTH==12)
saveRDS(Dec, "Intermediate data files/global_mean_abundance_from_ebird/december.RDS")
