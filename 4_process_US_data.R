library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidytable)
# library(epimdr2)
library(cowplot)
library(dplyr)
library(beepr)
library(fields)
library(rootSolve)
library(patchwork)
library(stringr)
library(readxl)

#### PROCESS US DATA -----------------------------------------------------------

# get FIPS codes by county
FIPS = read.csv("data/US-data/State__County_and_City_FIPS_Reference_Table_20251007.csv") %>%
  mutate(state_name = str_to_title(State.Name), 
         county_name = str_to_title(County.Name), 
         location_id = StCnty.FIPS.Code, 
         location_name = paste0(county_name, ", ", state_name)
  ) %>%
  select(state_name, county_name, location_id)

# get measles case
measles_cases = read.csv("data/US-data/measles_county_outbreaks_2026.02.13.csv") %>%
  mutate(location_name = paste0(county, ", ", state)) %>% 
  rename(county_name = county, 
         state_name = state) %>%
  left_join(unique(FIPS))

# do some manual updating of the data to associate each with a county
# get county-level detail on SC outbreak form dph.sc.gov (pulled on Feb 16, 2026)
# cherokee county is reported as <5 cases, so use 2 here (as floor(median))
# https://dph.sc.gov/diseases-conditions/infectious-diseases/measles-rubeola/measles-dashboard
upstate_sc = data.frame(county_name = c("Spartanburg", "Greenville", "Anderson", "Cherokee"), 
                        state_name = rep("South Carolina", 4), 
                        count = c(904, 35, 6, 2), 
                        location_id = c(45083, 45045, 45007, 45021)) %>%
  mutate(location_name = paste0(county_name, ", ", state_name))

# now add to measles cases
measles_cases = measles_cases %>% 
  filter(location_name != "Upstate, South Carolina") %>%
  bind_rows(upstate_sc)

# update single CT case to Fairfield county following https://portal.ct.gov/dph/home/newsroom/press-releases---2025/measles-case-in-connecticut?language=en_US
measles_cases = measles_cases %>%
  mutate(location_name = ifelse(location_name == "South Central Connecticut, Connecticut", "Fairfield, Connecticut", location_name)) %>%
  mutate(location_id = ifelse(location_name == "Fairfield, Connecticut", 9001, location_id)) %>% 
  mutate(county_name = ifelse(location_name == "Fairfield, Connecticut", "Fairfield", county_name))

# get measles vax
measles_vax = read.csv("data/US-data/mmr_data_us_counties.csv") %>%
  rename(location_id = FIPS,
         state_name = State,
         county_name = County) %>%
  mutate(location_name = paste0(county_name, ", ", state_name)) %>%
  mutate(mean_vax = mean(c(as.double(SY2017_18), as.double(SY2018_19), as.double(SY2019_20), as.double(SY2020_21), as.double(SY2021_22), as.double(SY2022_23), as.double(SY2023_24)), na.rm = TRUE), .by = c("location_name", "location_id", "county_name"))

# get birth rates and population sizes (2019, social explorer)
birth_pop = read.csv("data/US-data/social_explorer_2019_US_census_birth_pop.csv", stringsAsFactors = FALSE)[-1,] %>%
  mutate(county_name = gsub(" County", "", Name.of.Area), 
         state_name = substr(Qualifying.Name, gregexpr(", ", Qualifying.Name)[[1]][1] + 2, nchar(Qualifying.Name)), 
         .by = c("Qualifying.Name")) %>% 
  rename(total_pop = `Total.Population`, 
         births = `Births`, 
         birth_rate = `Births.Rate.per.1.000.Population`, 
         location_id = FIPS) %>%
  select(county_name, state_name, location_id, total_pop, births, birth_rate) %>%
  mutate(total_pop = as.integer(total_pop), 
         births = as.integer(births), 
         birth_rate = as.integer(birth_rate)) %>% 
  mutate(county_name = ifelse(county_name == "Do<f1>a Ana", "Dona Ana", county_name))

#### COMBINE DATASETS AND ADJUST MANUALLY WHERE NEEDED --------------------------------
us_data_corrected = measles_cases %>%
  summarize(total_cases = sum(count), .by = c("location_id", "location_name", "state_name")) %>%
  # some corrections for merging 
  mutate(location_name = ifelse(state_name == "District of Columbia", "District of Columbia, District of Columbia", location_name)) %>%
  mutate(location_id = ifelse(state_name == "District of Columbia", 11001, 
                              ifelse(location_name == "St. Johns, Florida", 12109, 
                                     ifelse(location_name == "McLennan, Texas", 48309, location_id)))) %>%
  full_join(
    birth_pop %>% 
      mutate(location_id = as.integer(location_id), 
             location_name = paste0(county_name, ", ", state_name)) %>%
      select(location_id, location_name, state_name, total_pop, births, birth_rate)
  ) %>%
  full_join(measles_vax %>% select(location_id, mean_vax)) %>%
  mutate(total_cases = ifelse(is.na(total_cases), 0, total_cases), 
         total_cases_per_pop = total_cases/total_pop) %>%
  mutate(log_mean_birth_rate = log(births/total_pop))# %>%
#filter(!is.na(location_id))

# now address states where cases were reported regionally rather than by county, so we cannot disaggregate
# when we can find the counties that belong to a given health region, we use a population-weighted average
# of birth rates and vax rates

# utah - map from https://files.epi.utah.gov/Utah%20measles%20dashboard.html then converted using UT county map
utah_conversion = data.frame(data_name = paste0(c(rep("Southwest Health District", 5), 
                                                  rep("Southeast Health District", 3), 
                                                  rep("Central", 6), 
                                                  #rep("TriCounty", 3),  exclude because no cases reported in this region
                                                  rep("Bear River", 3), 
                                                  rep("Weber-Morgan", 2)
), ", Utah"),
county_name = paste0(c("Beaver", "Iron", "Garfield", "Kane", "Washington", 
                       "Grand", "Emery", "Carbon", 
                       "Juab", "Millard", "Sanpete", "Sevier", "Piute", "Wayne", 
                       #"Daggett", "Uintah", "Duchesne", 
                       "Box Elder", "Cache", "Rich", 
                       "Weber", "Morgan"
), ", Utah"))


louisiana_conversion = data.frame(
  county_name = paste0(c("Livingston", "Tangipahoa", "St. Helena", "St. Tammany", "Washington", "Jefferson", "St. Bernard", "Orleans", "Plaquemines"), " Parish, Louisiana"),
  data_name = paste0(c(rep("Region 9", 5), rep("Region 1", 4)), ", Louisiana")
)


NYC_conversion = data.frame(
  county_name = paste0(c("Bronx", "Kings", "New York", "Queens", "Richmond"), ", New York"),
  data_name = "New York City, New York"
)

tennessee_conversion = data.frame(
  county_name = paste0(c("Cheatham", "Dickson", "Houston", "Humphreys", "Montgomery","Robertson", "Rutherford", "Stewart", "Sumner", "Trousdale", "Williamson", "Wilson", 
                         "Davidson", 
                         "Cannon", "Clay", "Cumberland", "DeKalb", "Fentress", "Jackson", "Macon", "Overton", "Pickett", "Putnam", "Smith", "Van Buren", "Warren"), ", Tennessee"),
  data_name = paste0(c(rep("Mid-Cumberland Region", 12), "Nashville-Davidson County Region", rep("Upper Cumberland Region", 13)), ", Tennessee")
)

virginia_conversion = read.csv("data/US-data/virginia_locality_hd.csv") %>%
  rename(county_name = Locality, 
         location_id = CountyFIPS, 
         data_name = HPR) %>%
  mutate(county_name = paste0(county_name, ", Virginia"), 
         data_name = paste0(data_name, ", Virginia")) %>%
  select(county_name, data_name)

conversion_all = bind_rows(utah_conversion,
                           louisiana_conversion, 
                           NYC_conversion, 
                           tennessee_conversion, 
                           virginia_conversion)

conversion_rates = us_data_corrected %>% filter(location_name %in% conversion_all$county_name) %>%
  left_join(conversion_all, by = c("location_name" = "county_name")) %>%
  summarize(mean_vax = weighted.mean(mean_vax, total_pop, na.rm = TRUE),
            births = sum(births),
            total_pop = sum(total_pop),
            .by = "data_name") %>%
  rename(location_name = data_name) %>%
  left_join(measles_cases %>%
              summarize(total_cases = sum(count), .by = c("location_id", "location_name", "state_name"))) %>%
  mutate(log_mean_birth_rate = log(births/total_pop), 
         total_cases_per_pop = total_cases/total_pop)

us_data_corrected = us_data_corrected %>% 
  filter(!(location_name %in% conversion_all$county_name), 
         !(location_name %in% conversion_rates$location_name)) %>%
  bind_rows(conversion_rates)

#### CHECK MERGING ---------------------------------------------------------------
tst_mu = seq(-5.5, -3.5, length.out = 20)#seq(log(1/45), log(1/200), length.out = 20)

n_counties_births <- birth_pop %>%
  summarize(n_counties_births = n_distinct(county_name), .by = "state_name")

n_counties_vax <- measles_vax %>%
  mutate(state_name = ifelse(state_name == "DC",
                             "District of Columbia",
                             state.name[match(state_name, state.abb)])) %>%
  summarize(n_counties_vax = n_distinct(county_name), .by = "state_name")

n_counties_cases <- measles_cases %>%
  summarize(n_counties_cases = n_distinct(county_name), .by = "state_name")

total_cases <- measles_cases %>%
  summarize(total_cases = sum(count), .by = "state_name")

total_cases_excluded <- us_data_corrected %>% 
  filter(is.na(mean_vax) | 
           log_mean_birth_rate < log(min(exp(tst_mu)))) %>%
  summarize(total_cases_excluded  = sum(total_cases), .by = "state_name")

# Merge the three summaries by state_name for comparison
county_counts_compare <- n_counties_births %>%
  full_join(n_counties_vax, by = "state_name") %>%
  full_join(n_counties_cases, by = "state_name") %>% 
  full_join(total_cases, by = "state_name") %>%
  full_join(total_cases_excluded, by = "state_name") %>%
  filter(state_name != "Puerto Rico")

# View(county_counts_compare)

# first look at how many cases occurred in states without vaccination data
states_no_vax = county_counts_compare %>% filter(is.na(n_counties_vax)) %>% pull(state_name)

#### KEEP DATA.FRAME OF NOTES --------------------------------------------------
# keep notes on states with mismatch between cases and vax data
notes_df = data.frame(
  state_name = c(
    "Alabama",
    "California",
    "Connecticut",
    "Iowa",
    "Kansas", 
    "Kentucky", 
    "Louisiana",
    "Minnesota",
    "New Jersey", 
    "New Mexico", 
    "New York", 
    "Oklahoma",
    "Oregon", 
    "Pennsylvania",
    "South Carolina",
    "Tennessee", 
    "Texas",
    "Utah",
    "Vermont",
    "Virginia", 
    states_no_vax
  ), 
  note = c(
    "1 case reported for North Alabama, could not identify county", 
    "4 cases in unknown county", 
    "1 case reported in South Central Connecticut which was identified as Fairfield county following https://portal.ct.gov/dph/home/newsroom/press-releases---2025/measles-case-in-connecticut?language=en_US", 
    "9 cases reported in Central Eastern and Western Iowa, but no reliable county decomposition available", 
    "20 cases reported in unknown county", 
    "1 case reported in unknown county",
    "health dept regions manually mapped to counties from https://ldh.la.gov/oph-regional-offices",
    "10 cases reported in Twin Cities Metro Area (and not clear which counties), 12 cases reported in unknown county", 
    "3 cases reported in unknown county",
    "6 cases reported in Sandoval, but no vax data available",
    "20 cases reported in New York City and population weighted average of vax and birth rates for NYC counties used (Bronx, Kings, New York, Queens, Richmond)",
    "18 cases (all cases) reported in unknown county",
    "2 cases reported in unknown county",
    "1 case reported in unknown county",
    "retrieved county-level detail on from dph.sc.gov on Feb 16, 2026; https://dph.sc.gov/diseases-conditions/infectious-diseases/measles-rubeola/measles-dashboard",
    "cases reported in public health regions, manually converted with https://www.tn.gov/health/health-program-areas/oralhealth/clinics/oral-health-regions-contacts.html",
    "1 case reported in unknown county",
    "cases reported by region so aggregated using https://files.epi.utah.gov/Utah%20measles%20dashboard.html then converted using UT county map",
    "1 case reported in unknown county",
    "cases reported regionally and converted from https://www.vdh.virginia.gov/content/uploads/2020/05/Locality-to-HD-to-HPR.pdf", 
    rep("no vaccination data available", length(states_no_vax))
  )
)

#### SUMMARIZE INTO SINGLE TABLE (TABLE S1) ------------------------------------
# okay but this looks at each independently, so also want to check
# that we have the same counties in each dataset
completeness_check = us_data_corrected %>% 
  mutate(cases_flag = total_cases > 0, 
         vax_flag = !is.na(mean_vax), 
         case_vax_flag_comb = ifelse(cases_flag & vax_flag, "case reported and vax data", 
                                     ifelse(cases_flag & !vax_flag, "case reported, no vax data", 
                                            ifelse(!cases_flag & vax_flag, "no case reported", "no case reported")))) %>%
  left_join(birth_pop %>% 
              mutate(location_name = paste0(county_name, ", ", state_name)) %>%
              select(state_name, location_name)) %>%
  summarize(n = n_distinct(location_name), .by = c("case_vax_flag_comb", "state_name")) %>%
  pivot_wider(names_from = case_vax_flag_comb, values_from = n) %>%
  left_join(county_counts_compare) %>%
  filter(state_name != "Puerto Rico") %>% 
  mutate(pct_cases_excluded = total_cases_excluded/total_cases) %>%
  left_join(notes_df)
write.csv(completeness_check, "data/output-data/county_data_completeness_check.csv", row.names = FALSE)

# save object
saveRDS(us_data_corrected, "data/output-data/us_data_corrected.rda")

