### In this script, we download shapefiles for North America and for bodies of water,
### which we'll use later for generating a habitat data for each eBird checklist and 
### for making maps.

### As with the other scripts in this folder, much of this code has been taken from 
### the "Best Practices for Using eBird Data" guide by Strimas-Mackey et al. (2020):

### Strimas-Mackey, M., W.M. Hochachka, V. Ruiz-Gutierrez, O.J. Robinson, E.T. Miller,
### T. Auer, S. Kelling, D. Fink, A. Johnston. 2020. Best Practices for Using eBird 
### Data. Version 1.0. https://cornelllabofornithology.github.io/ebird-best-practices/. 
### Cornell Lab of Ornithology, Ithaca, New York. https://doi.org/10.5281/zenodo.3620739

library(tidyverse)
library(sf)
library(rnaturalearth)

setwd("/Volumes/project/CarlingLab/pdoughe1/Passerina_hybrid_zone_changes_parental_abudnance_estimation") # accessing eBird data in Alcova

#=============================================================================================
#     downloading shapefiles
#=============================================================================================

# file to save shapefiles
gpkg_dir <- "spatial_data"
if (!dir.exists(gpkg_dir)) {
  dir.create(gpkg_dir)
}
f_ne <- file.path(gpkg_dir, "gis_data.gpkg")

# downloading shapefile of North America
ne_land <- ne_download(scale = 50, category = "cultural",
                       type = "admin_0_countries_lakes",
                       returnclass = "sf") %>%
  filter(CONTINENT == "North America") %>%
  st_set_precision(1e6) %>%
  st_union()
write_sf(ne_land, f_ne, "ne_land")

# also downloading shapefiles for rivers and lakes, which we'll add to our maps for reference
ne_rivers <- ne_download(scale = 50, type = 'rivers_lake_centerlines', 
                         category = 'physical',
                         returnclass = "sf") %>%
  st_set_precision(1e6) %>%
  st_union()
write_sf(ne_rivers, f_ne, "ne_rivers")

ne_lakes <- ne_download(scale = 50, type = 'lakes', 
                        category = 'physical',
                        returnclass = "sf") %>%
  filter(adm0_a3 %in% c("USA", "CAN", "MEX")) %>%
  st_set_precision(1e6) %>%
  st_union()
write_sf(ne_lakes, f_ne, "ne_lakes")