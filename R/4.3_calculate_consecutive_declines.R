library(ggplot2)
library(dplyr)
library(reshape2)
library(readxl)
library(cowplot)

#### FUNCTION ------------------------------------------------------------------
get_max_drop = function(vacc, data_start_yr, data_end_yr = 2024){
  n_yrs = length(vacc[!is.na(vacc)])
  # vacc_interp = interp_vals(vacc)
  vacc_interp = ifelse(vacc > 1, 1, vacc)
  if(n_yrs <= 1){
    return(data.frame(nyears_data = n_yrs,
                      drop = NA, 
                      nyears_drop = NA, 
                      start_yr = NA))
  }
  d = diff(vacc_interp)
  drops = rep(NA, length(d))
  drop_size = rep(NA, length(d))
  for(i in 1:length(d)){
    if(d[i] < 0 & !is.na(d[i])){
      counter = 1
      drops[i] = d[i]
      drop_size[i] = 1
      while(d[i+counter]< 0 & !is.na(d[i+counter])){
        drops[i] = drops[i] + d[i+counter]
        drop_size[i] = drop_size[i] + 1
        counter = counter + 1
      }
    }
  }
  if(all(is.na(drops))){
    return(data.frame(nyears_data = n_yrs,
                      drop = NA, 
                      nyears_drop = NA, 
                      start_yr = NA))
  }
  else{
    biggest_drop = min(drops, na.rm = TRUE)
    which_biggest_drop = which.min(drops)
  }
  return(data.frame(nyears_data = n_yrs,
                    drop = biggest_drop, 
                    nyears_drop = drop_size[which_biggest_drop], 
                    start_yr = (data_start_yr:data_end_yr)[which.min(drops)]))
}

#### COVERAGE DROPS IN WHO DATA  -----------------------------------------------
folder = "data/WHO-data/"

start_yr = 2000

who_vacc = read.csv(paste0(folder, "measlesVaccCoverFirstDose.csv")) %>%
  filter(COVERAGE_CATEGORY == "WUENIC") %>%
  rename(year = YEAR, coverage = COVERAGE) %>% 
  mutate(year = as.integer(year)) %>%
  rename(country_name = NAME) %>% 
  select(country_name, year, coverage) %>% 
  mutate(coverage = coverage/100)

# exclude 3 countries where vax data started after 2000
who_vacc_sub = who_vacc %>% 
  filter(year >= start_yr, !(country_name %in% c("South Sudan", "Montenegro", "Timor-Leste")))

who_drops_by_country = who_vacc_sub %>% 
  arrange(country_name, year) %>%
  dplyr::reframe(get_max_drop(coverage, start_yr), .by = c("country_name"))

bin_width = 0.02
drop_bins = seq(-1, 0, bin_width)
who_drops_by_country_summ  = who_drops_by_country %>%
  mutate(row_id = seq_len(n())) %>%
  mutate(drop_bin = ifelse(!is.na(drop), drop_bins[min(which(drop < drop_bins))], NA), .by = c("row_id"))

saveRDS(who_drops_by_country_summ, "data/output-data/drops_by_country_WHO.rda")

#### COVERAGE DROPS IS US DATA -------------------------------------------------
# FIPS: https://data.transportation.gov/Railroads/State-County-and-City-FIPS-Reference-Table/eek5-pv8d/about_data

FIPS = read.csv("data/US-data/State__County_and_City_FIPS_Reference_Table_20251007.csv") %>%
  select(State.Name, State.Code, County.Name, StCnty.FIPS.Code) %>%
  rename(state_name = State.Name, state_abbrev = State.Code,
         county_name = County.Name, location_id = StCnty.FIPS.Code
  ) %>%
  mutate(state_name = tools::toTitleCase(tolower(state_name)), 
         county_name = tools::toTitleCase(tolower(county_name))) %>% 
  unique()

vacc = read.csv("data/US-data/mmr_data_us_counties.csv") %>% 
  melt(c("FIPS", "County", "State")) %>% 
  mutate(start_year = as.integer(substr(variable, 3,6)), 
         value = as.double(value)) %>%
  rename(location_id = FIPS, state_abbrev = State, county_name = County)

drops_by_county = vacc %>% 
  arrange(state_abbrev, location_id, start_year) %>%
  dplyr::reframe(get_max_drop(value, 2017, 2023), .by = c("county_name", "location_id", "state_abbrev"))

drops_by_county_summ  = drops_by_county %>%
  mutate(row_id = seq_len(n())) %>%
  mutate(drop_bin = ifelse(!is.na(drop), drop_bins[min(which(drop <= drop_bins))], NA), .by = c("row_id"))

saveRDS(drops_by_county_summ, "data/output-data/drops_by_county_US.rda")