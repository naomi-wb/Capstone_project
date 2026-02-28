##

library(geojsonR)
library(geojsonio)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(tibble)
library(sf)
library(ggplot2)

## util function so it can be applied to any grid shapefile ---- 
get_shapefile_res <- function(grid) {
  grid_cells <- group_by(grid, centroid_lng) %>% 
    arrange(centroid_lng, centroid_lat) 
  
  grid_cells[2 ,]$centroid_lat - grid_cells[1 ,]$centroid_lat
}

## data ---- 

grid_cells_sp <- geojson_read("Data/Spatial_data/grid_5_degree.geojson", what = "sp") 
grid_cells_sf <- as(grid_cells_sp, "sf")

world_sp <- ne_countries(returnclass = "sp")
world_sf <- ne_countries(returnclass = "sf")

## intersection, sp types are used to get projected area in addition to true area ---- 

intersection <- raster::intersect(grid_cells_sp, world_sp) 

## area as it appears on the grid, not true area ---- 
grid_cell_projected_area <- get_shapefile_res(grid_cells_sf) ^ 2 

intersection_data <- cbind(intersection@data, projected_area = unlist(lapply(intersection@polygons, function(x) {x@area}))) %>% 
  dplyr::select(ID, centroid_lat, centroid_lng, area_km2, projected_area) %>%
  group_by(ID) %>% 
  mutate(projected_area_sum = sum(projected_area), 
         terrestrial_land = projected_area_sum/grid_cell_projected_area * area_km2) %>% 
  dplyr::select(ID, terrestrial_land) %>% 
  unique()

grid_cells_sf <- left_join(grid_cells_sf, intersection_data, 
                           by = c("ID" = "ID"), 
                           all.x = TRUE) 
grid_cells_sf$terrestrial_land[is.na(grid_cells_sf$terrestrial_land)] <- 0
grid_cells_sf <- mutate(grid_cells_sf, percentage = (terrestrial_land/area_km2) * 100)

saveRDS(grid_cells_sf, "imputation_analysis/grid_percent_land.RDS")

ggplot() +
 geom_sf(data = grid_cells_sf, aes(fill = percentage)) +
 scale_colour_gradient(low = "blue", high = "green") +
 scale_fill_viridis_c(option = "magma") +
 geom_sf(data = world_sf, colour = "white", fill = NA)