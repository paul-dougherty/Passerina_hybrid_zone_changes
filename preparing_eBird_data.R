##testing 

### In this script, we use the auk package to filter eBird checklists 
### down to only those that we can use to predict the distributions/
### relative abundances of focal taxa during the breeding season. Our 
### filtering decisions were largely inspired by Strimas-Mackey et al.
### (2020), although we also considered the life history of focal taxa
### (see below). For the sake of clarity, I filtered checklists for each
### taxon separately (as opposed to using a for loop or something similar).

### As with the other scripts in this folder, much of this code has been taken from 
### the "Best Practices for Using eBird Data" guide by Strimas-Mackey et al. (2020):

### Strimas-Mackey, M., W.M. Hochachka, V. Ruiz-Gutierrez, O.J. Robinson, E.T. Miller,
### T. Auer, S. Kelling, D. Fink, A. Johnston. 2020. Best Practices for Using eBird 
### Data. Version 1.0. https://cornelllabofornithology.github.io/ebird-best-practices/. 
### Cornell Lab of Ornithology, Ithaca, New York. https://doi.org/10.5281/zenodo.3620739

#packages for preparing eBird data
library(tidyverse)
library(auk)
library(lubridate)

# defining the directory for the eBird data
ebd_dir <-"/Volumes/commons/CarlingLab/eBird Data/Data for looking at relative abundance"
#ebd_dir <- "/pfs/tsfs1/gscratch/pdoughe1" # if doing on Teton


#=============================================================================================
#     preparing eBird data
#=============================================================================================
## Download eBird data for focal taxa from https://ebird.org/science/use-ebird-data 
## (after obtaining permission, of course). While you can start with the entire
## eBird database and then filter down to taxa of interest, as we did when 
## estiamting overall hybridization rates, it's much easier/faster to download 
## only records for these taxa. In this script, I downloaded eBird records for 
## Indigo Bunting, Lazuli Bunting, and Indigo x Lazuli buntings.

# filtering checklists for first parental species -------------------------
# filtering checklists for Indigo Bunting

ebd <- auk_ebd("ebd_indbun_relMar-2020.txt",
               file_sampling = "ebd_sampling_relMar-2020.txt")

ebd_filters <- ebd %>% 
  auk_species(species = c("Indigo Bunting")) %>%
  ## restrict to the standard traveling and stationary count protocols
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  auk_date(date = c("*-06-01", "*-07-16")) %>% #June and July, but only up to July 16th as Lazuli Bunting departs breeding grounds in late July
  auk_last_edited(date = c("2010-01-01" , "2018-12-31")) %>% #2010 to 2018 just like he Justyn et al. paper
  auk_country(country = c("United States", "Canada", "Mexico")) %>%
  auk_complete()
ebd_filters # viewing the filters

## merging eBird and sampling data for one species at a time
f_ebd <- file.path(ebd_dir, "ebd_indigo_bunting.txt")
f_sampling <- file.path(ebd_dir, "ebd_indigo_bunting_sampling.txt")

## filtering the zero-filled data, only run if the files don't already exist,
auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling, overwrite = TRUE)

#reading in zero-filled data into R
if (!exists("indigo_bunting_zf")) {
  indigo_bunting_zf <- auk_zerofill(f_ebd, f_sampling, collapse = TRUE)
}

## cleaning up eBird dataset and preparing it for analyses
# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

## cleaning up variables
indigo_bunting_zf <- indigo_bunting_zf %>% 
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

# additional filtering
indigo_zf_filtered <- indigo_bunting_zf %>% 
  filter(
    # effort filters
    duration_minutes <= 5 * 60,
    effort_distance_km <= 5,
    # last 10 years of data
    year >= 2010,
    # 10 or fewer observers
    number_observers <= 10)

# filtering down to relevant variables
indigo_bunting_pred <- indigo_zf_filtered %>% 
  select(checklist_id, observer_id, sampling_event_identifier,
         scientific_name,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers)

## saving to csv file
write.csv(indigo_bunting_pred, "indigo_bunting_pred_mx_half_july.csv", na = "", row.names=FALSE)


# filtering checklists for second parental species -------------------------
# filtering checklists for Lazuli Bunting
ebd <- auk_ebd("ebd_lazbun_relMar-2020.txt", 
               file_sampling = "ebd_sampling_relMar-2020.txt")

ebd_filters <- ebd %>% 
  auk_species(species = c("Lazuli Bunting")) %>%
  ## restrict to the standard traveling and stationary count protocols
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  auk_date(date = c("*-06-01", "*-07-16")) %>% #June and July, only up to July 16th
  auk_last_edited(date = c("2010-01-01" , "2018-12-31")) %>%
  auk_country(country = c("United States", "Canada", "Mexico")) %>%
  auk_complete()
ebd_filters # viewing the filters

## merging eBird and sampling data for one species at a time
f_ebd <- file.path(ebd_dir, "ebd_lazuli_bunting.txt")
f_sampling <- file.path(ebd_dir, "ebd_lazuli_bunting_sampling.txt")

## filtering the zero-filled data, only run if the files don't already exist,
auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling, overwrite = TRUE)

#reading in zero-filled data into R
if (!exists("lazuli_bunting_zf")) {
  lazuli_bunting_zf <- auk_zerofill(f_ebd, f_sampling, collapse = TRUE)
}

## cleaning up eBird dataset and preparing it for analyses

## cleaning up variables
lazuli_bunting_zf <- lazuli_bunting_zf %>% 
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

# additional filtering
lazuli_zf_filtered <- lazuli_bunting_zf %>% 
  filter(
    # effort filters
    duration_minutes <= 5 * 60,
    effort_distance_km <= 5,
    # last 10 years of data
    year >= 2010,
    # 10 or fewer observers
    number_observers <= 10)

# filtering down to relevant variables
lazuli_bunting_pred <- lazuli_zf_filtered %>% 
  select(checklist_id, observer_id, sampling_event_identifier,
         scientific_name,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers)

## saving to csv file
write.csv(lazuli_bunting_pred, "lazuli_bunting_pred_mx_half_july.csv", na = "", row.names=FALSE)