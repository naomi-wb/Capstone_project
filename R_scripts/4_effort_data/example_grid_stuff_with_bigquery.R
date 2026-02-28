# summarize the total number of minutes spent sampling in each grid
# all done in bigquery

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


# create connection with online database
con <- DBI::dbConnect(bigrquery::bigquery(),
                      dataset= "ebird",
                      project="ebird-database",
                      billing="ebird-database")

# create ebird table
ebird <- tbl(con, 'ebird_qa')

# create locality_id and grid_id table
grids <- tbl(con, 'grid_5_joined_with_sites')


dat <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, LOCALITY_ID, COMMON_NAME, 
                OBSERVATION_COUNT, OBSERVATION_DATE, SCIENTIFIC_NAME, DURATION_MINUTES) %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  group_by(SAMPLING_EVENT_IDENTIFIER, LOCALITY_ID) %>%
  summarize(checklist_minutes=sum(DURATION_MINUTES)) %>%
  left_join(., grids, by="LOCALITY_ID") %>%
  group_by(grid_id) %>%
  summarize(total_minutes=sum(checklist_minutes)) %>%
  collect(n=Inf)


# download some natural earth data
countries <- ne_countries(returnclass = "sf") %>%
  st_geometry()

graticules <- ne_download(type = "graticules_15", category = "physical",
                          returnclass = "sf") %>%
  st_geometry()

bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>%
  st_geometry()

# read in underlying map of grids for 5 degree grid cells
grid_5_degree <- st_read("Data/Spatial_data/grid_5_degree.geojson")


# combine dat from above which is a summary from bigquery
# with the 5 degree grids
grid_minutes_summary <- grid_5_degree %>%
  rename(grid_id=ID) %>%
  left_join(., dat, by="grid_id") %>%
  replace_na(list(total_minutes=0)) %>%
  mutate(minutes_by_area=total_minutes/area_km2)


# make a plot of the grids over the world
plasma_pal <- c(viridis::plasma(n = 10))

ggplot()+
  geom_sf(data = bb, col = "grey20", fill = "transparent")+
  geom_sf(data = countries, fill = "black", col = "grey40", lwd = 0.3)+
  geom_sf(data= grid_minutes_summary, aes(fill=log(total_minutes)), alpha=0.8)+
  theme_minimal()+
  scale_fill_gradientn(colours=plasma_pal)+
  theme(axis.text = element_blank())+
  theme(legend.position="bottom")+
  ggtitle("Total number of minutes sampling per grid cell in the world")


# read in clements data
clements <- read_csv("Data/Clements-Checklist-v2018-August-2018.csv") %>%
  dplyr::filter(category == "species") %>%
  rename(COMMON_NAME = `English name`) %>%
  rename(SCIENTIFIC_NAME = `scientific name`) %>%
  dplyr::select(COMMON_NAME, SCIENTIFIC_NAME, order, family)
