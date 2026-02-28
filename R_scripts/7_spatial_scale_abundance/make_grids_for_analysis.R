# This is an R script to make grids of the world
# at different degrees (1, 5, 10, 20)
# and then save these grids out to the repository as a geojson
# then in order to use these geojson polygons in BigQuery
# they must be uploaded as a csv file
# and this is done in the command line using
# the code at the bottom of this script

# packages
library(sf)
library(rnaturalearth)
library(mapproj)
library(ggplot2)
library(dplyr)

# download some natural earth data
countries <- ne_countries(returnclass = "sf") %>%
  st_geometry()

graticules <- ne_download(type = "graticules_15", category = "physical",
                          returnclass = "sf") %>%
  st_geometry()

bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>%
  st_geometry()

# make a simple plot of the world
# will use this for later
ggplot()+
  geom_sf(data = bb, col = "grey20", fill = "transparent")+
  geom_sf(data = graticules, col = "grey20", lwd = 0.1)+
  geom_sf(data = countries, fill = "grey80", col = "grey40", lwd = 0.3)+
  theme_minimal()+
  theme(axis.text = element_blank())

# now make a 5x5 degree grid covering the world
globe_grid_5 <- st_make_grid(bb, cellsize = c(5, 5),
                                  crs = 4326, what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('ID' = 1:length(.)))

globe_grid_5 <- globe_grid_5 %>%
  mutate(area_km2=as.numeric(st_area(.))/1000) %>%
  bind_cols(., as.data.frame(st_coordinates(st_centroid(.)))) %>%
  rename(centroid_lng=Y) %>%
  rename(centroid_lat=X)

# now make a 10x10 degree grid covering the world
globe_grid_10 <- st_make_grid(bb, cellsize = c(10, 10),
                             crs = 4326, what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('ID' = 1:length(.)))

globe_grid_10 <- globe_grid_10 %>%
  mutate(area_km2=as.numeric(st_area(.))/1000) %>%
  bind_cols(., as.data.frame(st_coordinates(st_centroid(.)))) %>%
  rename(centroid_lng=Y) %>%
  rename(centroid_lat=X)

# now make a 15x15 degree grid covering the world
globe_grid_15 <- st_make_grid(bb, cellsize = c(15, 15),
                             crs = 4326, what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('ID' = 1:length(.)))

globe_grid_15 <- globe_grid_15 %>%
  mutate(area_km2=as.numeric(st_area(.))/1000) %>%
  bind_cols(., as.data.frame(st_coordinates(st_centroid(.)))) %>%
  rename(centroid_lng=Y) %>%
  rename(centroid_lat=X)

# now make a 20x20 degree grid covering the world
globe_grid_20 <- st_make_grid(bb, cellsize = c(20, 20),
                             crs = 4326, what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('ID' = 1:length(.)))

globe_grid_20 <- globe_grid_20 %>%
  mutate(area_km2=as.numeric(st_area(.))/1000) %>%
  bind_cols(., as.data.frame(st_coordinates(st_centroid(.)))) %>%
  rename(centroid_lng=Y) %>%
  rename(centroid_lat=X)

# write out each of these as a geojson file
st_write(globe_grid_5, "Data/Spatial_data/grid_5_degree.geojson")
st_write(globe_grid_10, "Data/Spatial_data/grid_10_degree.geojson")
st_write(globe_grid_15, "Data/Spatial_data/grid_15_degree.geojson")
st_write(globe_grid_20, "Data/Spatial_data/grid_20_degree.geojson")


# an example of a plot with the 5 degree grids drawn
ggplot()+
  geom_sf(data = bb, col = "grey20", fill = "transparent")+
  geom_sf(data = graticules, col = "grey20", lwd = 0.1)+
  geom_sf(data = countries, fill = "grey80", col = "grey40", lwd = 0.3)+
  geom_sf(data=globe_grid_5, aes(fill="green"), alpha=0.7)+
  theme_minimal()+
  theme(axis.text = element_blank())


# These are commands to copy into the command prompt
# which used ogr2ogr from GDAL to
# convert from the geojson to a csv
# where the geography is a separate column


#ogr2ogr -f csv -dialect sqlite -sql "select AsGeoJSON(geometry) AS geom, * from grid_5_degree" C:\Users\CTC\Desktop\Research\Current_projects\how_many_birds\Data\Spatial_data\grid_5_degree.csv C:\Users\CTC\Desktop\Research\Current_projects\how_many_birds\Data\Spatial_data\grid_5_degree.geojson
#ogr2ogr -f csv -dialect sqlite -sql "select AsGeoJSON(geometry) AS geom, * from grid_10_degree" C:\Users\CTC\Desktop\Research\Current_projects\how_many_birds\Data\Spatial_data\grid_10_degree.csv C:\Users\CTC\Desktop\Research\Current_projects\how_many_birds\Data\Spatial_data\grid_10_degree.geojson
#ogr2ogr -f csv -dialect sqlite -sql "select AsGeoJSON(geometry) AS geom, * from grid_15_degree" C:\Users\CTC\Desktop\Research\Current_projects\how_many_birds\Data\Spatial_data\grid_15_degree.csv C:\Users\CTC\Desktop\Research\Current_projects\how_many_birds\Data\Spatial_data\grid_15_degree.geojson
#ogr2ogr -f csv -dialect sqlite -sql "select AsGeoJSON(geometry) AS geom, * from grid_20_degree" C:\Users\CTC\Desktop\Research\Current_projects\how_many_birds\Data\Spatial_data\grid_20_degree.csv C:\Users\CTC\Desktop\Research\Current_projects\how_many_birds\Data\Spatial_data\grid_20_degree.geojson





