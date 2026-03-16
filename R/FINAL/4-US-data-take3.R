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

source("R/FINAL/0_helper-functions.R")
source("R/FINAL/0_SIR-age-functions.R")

#### MODEL SETUP ---------------------------------------------------------------
source("R/FINAL/2.1_setup-WAIFW.R") # this will add waifw to environment, which is list of WAIFW matrices to test
n_waifw = length(waifw)

compartments = c("S", "I", "R")

R0 <- 17
paras = c(mu = 1/80, N = 500000, beta1 = 0,
          gamma = 365/14, delta = 0, p = 0, phi = 1) # delta = 1e-4
paras["beta0"]  = (paras["gamma"] + paras["mu"])*R0

tst_mu = seq(-5.5, -3.5, length.out = 20)#seq(log(1/45), log(1/200), length.out = 20)

start_vax = 0.95

#### GET R0 VALUES FOR EACH WAIFW ----------------------------------------------
find_scalar = function(s, R0, waifw, S, beta0, gamma, mu, N){
  # print(paste0("s: ", s, " R0: ", get_Rt(waifw, S, beta0*s, gamma, N)))
  diff = get_Rt(waifw, S, beta0*s, gamma, mu, N) - R0
  return(abs(diff))
}

scalars_all = vector("list", length(tst_mu))
vax_equilib_long = vector("list", length(tst_mu))
for(j in 1:length(tst_mu)){
  tmp_mu = exp(tst_mu[j])
  fert = rep(tmp_mu, length(age_classes))
  mort = rep(tmp_mu, length(age_classes))
  stable_age = findStableStruct(age_classes, mort, fert, 1/52)$stable.age
  scalars_tmp = data.frame(waifw_id = 5, scalar = NA, diff = NA)
  i = 5
  # for(i in 5:length(waifw)){
    o = optimize(f = find_scalar, tol = 1e-8, interval = c(0, 100), R0 = R0,
                 waifw = waifw[[5]], S = stable_age, beta0 = paras["beta0"], gamma = paras["gamma"],  mu = tmp_mu, N = 1)
    print(get_Rt(waifw[[5]], stable_age, paras["beta0"]*o$minimum, paras["gamma"],  mu = tmp_mu, N = 1))
    scalars_tmp[1, 2:3] = c(o$minimum, o$objective)
    # scalars_tmp[i, 2:3] = c(o$minimum, o$objective)
  # }
  scalars_all[[j]]  = scalars_tmp
  # pre-vax equilibrium = (1-vax_cov)*N(a) where N(a) is stable age distribution
  vax_equilib_long[[j]] = expand.grid(age = age_classes, 
                                 start_vax = start_vax, 
                                 variable = c("S", "I", "R")) %>%
    left_join(data.frame(age = age_classes, 
                         bin_width = bin_width, 
                         N = stable_age*paras["N"])) %>%
    mutate(value = ifelse(variable == "S", N*(1-start_vax), ifelse(variable == "I", 0, N*start_vax)))
}

# create a list of parameters for 
paras_all = lapply(1:length(tst_mu), function(i){paras_tmp = paras; paras_tmp["beta0"] = paras_tmp["beta0"]*scalars_all[[i]][,2]; paras_tmp["mu"] = exp(tst_mu[i]); return(paras_tmp)})

#### CALCULATE HONEYMOON TIME --------------------------------------------------
# i.e., simulate across the full range of release vax values
release_vax_full = seq(0.5, 1, 0.02)
nyears_postrelease = 20
chosen_dt = 1/52

# generate data.frame of simulations to run
# note: do not need to repeat for all WAIFWs because there are no new infections
susc_after_release_full = vector("list", length(tst_mu))
for(j in 1:length(tst_mu)){
  print(paste0("j: ", j, "/", length(tst_mu)))
  tmp_mu = exp(tst_mu[j])
  fert = rep(tmp_mu, length(age_classes))
  mort = rep(tmp_mu, length(age_classes))
  tst_release_full = expand.grid(start_vax = start_vax, 
                                 release_vax = release_vax_full)
  tmp = vector("list", nrow(tst_release_full))
  for(i in 1:nrow(tst_release_full)){
    tmp_start_vax = tst_release_full[i, "start_vax"]
    tmp_release_vax = tst_release_full[i, "release_vax"]
    IC_manual = vax_equilib_long[[j]] %>% filter(start_vax == tmp_start_vax) %>% arrange(age)
    names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
    IC_manual = IC_manual$value
    names(IC_manual) = names_IC
    tmp[[i]] = run_ode(
      age_classes = age_classes, mort = mort, fert = fert, start_pop = paras["N"],
      compartments = compartments, dt = chosen_dt, waifw = waifw[[5]],
      vax_change_times = c(0), vax_rates = c(tmp_release_vax), 
      IC_type = "manual", IC_manual = IC_manual, max_t = nyears_postrelease, 
      params = paras_all[[j]],
      adjust_beta_flag = FALSE, plot_flag = FALSE
    )
  }
  susc_after_release_full[[j]] = bind_rows(tmp, .id = "sim_id") %>%
    mutate(sim_id = as.integer(sim_id)) %>%
    left_join(tst_release_full %>% mutate(sim_id = seq_len(dplyr::n())))
}
beep()

susc_after_release_long_full = bind_rows(susc_after_release_full, .id = "mu_id") %>%
  mutate(mu_id = as.integer(mu_id)) %>%
  left_join(bind_rows(lapply(paras_all, function(i){data.frame(t(i))}) ) %>% mutate(mu_id = seq_len(length(tst_mu)))) 

# POLYMOD ONLY
rt_after_release_full = susc_after_release_long_full %>% 
  filter(variable == "S") %>%
  summarize(Rt = get_Rt(waifw = waifw[[5]], S = value, beta0 = beta0, gamma = gamma, N = N, mu = mu), 
            .by = c("time", "start_vax", "release_vax", "mu_id", "mu"))

honeymoon_period = rt_after_release_full %>%
  left_join(bind_rows(lapply(paras_all, function(i){data.frame(t(i))}) ) %>% mutate(mu_id = seq_len(length(tst_mu)))) %>%
  filter(Rt > 1) %>%
  mutate(min_time = min(time), .by = c("mu_id", "start_vax", "release_vax", "mu")) %>%
  filter(time == min_time) %>%
  select(-min_time)

# saveRDS(honeymoon_period, "R/FINAL/data/honeymoon_period_by_birthrate.rds")

#### NOW GET REAL DATA ---------------------------------------------------------

# get FIPS codes by county
FIPS = read.csv("data/US-data/State__County_and_City_FIPS_Reference_Table_20251007.csv") %>%
  mutate(state_name = str_to_title(State.Name), 
         county_name = str_to_title(County.Name), 
         location_id = StCnty.FIPS.Code, 
         location_name = paste0(county_name, ", ", state_name)
  ) %>%
  select(state_name, county_name, location_id)

# get measles case
measles_cases_old = read.csv("data/US-data/measles_county_all_updates.csv")
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
plt_df = measles_cases %>%
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

conversion_rates = plt_df %>% filter(location_name %in% conversion_all$county_name) %>%
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

plt_df = plt_df %>% 
filter(!(location_name %in% conversion_all$county_name), 
!(location_name %in% conversion_rates$location_name)) %>%
  bind_rows(conversion_rates)



#### PLOT --------------------------------------------------------------------------


plot_honeymoon = expand.grid(mu = exp(tst_mu), 
                             release_vax = release_vax_full) %>%
  left_join(honeymoon_period %>% select(mu, release_vax, time))

plt_df %>% 
  filter(births > 0) %>%
  pull(log_mean_birth_rate) %>%
  range(na.rm = TRUE)

us_dat_fig = ggplot(plt_df %>% filter(log_mean_birth_rate > log(min(plot_honeymoon$mu)))) + 
  geom_tile(data = plot_honeymoon, aes(x = log(mu), y = release_vax, fill = time)) +
  # geom_contour(data = plot_honeymoon, aes(x = log(mu), y = release_vax, z = time),
  #              color = "black", breaks = c(1, 3, 5, 7), linewidth = 0.2) +
  # metR::geom_text_contour(data = plot_honeymoon,
  #                          aes(x = log(mu), y = release_vax, z = time),
  #                          breaks = c(1, 3, 5, 7), size = 3) +
  geom_point(aes(x = log_mean_birth_rate, y = as.double(mean_vax)), alpha = 0.8, shape = 21, color = "darkgray", size = 0.5, stroke = 0.2) +
  geom_point(data = plt_df %>% filter(total_cases_per_pop > 0, log_mean_birth_rate > log(min(plot_honeymoon$mu))),
             aes(x = log_mean_birth_rate, y = as.double(mean_vax), size = total_cases_per_pop, color = log(total_cases_per_pop*1e4))) +
  ggrepel::geom_text_repel(data = plt_df %>% filter(total_cases_per_pop > 0, log_mean_birth_rate > log(min(plot_honeymoon$mu))) %>% 
                             filter(total_cases_per_pop > quantile(total_cases_per_pop, 0.95)),
                           aes(x = log_mean_birth_rate, y = as.double(mean_vax), label = location_name),
                           box.padding = unit(1.1, "lines"),
                           point.padding = unit(0.7, "lines"),
                           min.segment.length = 0,
                           segment.size = 0.4,
                           segment.color = 'black', size = 2
  ) +
  guides(size = "none") +
  scale_color_viridis_c(breaks = log(c(1, 10, 100)), #c(-2, 0, 2, 4, 6), 
                       labels = c(1, 10, 100),
                       #trans = scales::pseudo_log_trans(sigma = 0.001), 
                       na.value = "darkgray", name = "cases per\n10,000") +
  scale_fill_viridis_c(option = "rocket", na.value = "#FAEBDDFF", name = "time to\nRe > 1") +
  scale_size_continuous(range = c(0,3)) +
  scale_x_continuous(expand = c(0,0),
                     breaks = c(log(c(50, 100, 200)/1e4)), 
                     labels = c(50, 100, 200),
                     name = "births per 10,000") +
  scale_y_continuous(expand = c(0,0), name = "vaccination coverage", 
                     labels = scales::percent) +
  theme_bw(base_size = 7) + 
  theme(legend.key.width = unit(0.5, "cm"),
    legend.position = "bottom")

save("us_dat_fig",
     file="R/FINAL/data/us_dat_fig.rda")

#### CHECK MERGING ---------------------------------------------------------------
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

total_cases_excluded <- plt_df %>% 
filter(is.na(mean_vax) | 
       log_mean_birth_rate < log(min(plot_honeymoon$mu))) %>%
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

# plt_df %>% filter(state_name %in% states_no_vax) %>%
  # summarize(total_cases = sum(total_cases), .by = "state_name") %>%
  # View()

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

# okay but this looks at each independently, so also want to check
# that we have the same counties in each dataset
completeness_check = plt_df %>% 
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
write.csv(completeness_check, "R/FINAL/data/county_data_completeness_check.csv", row.names = FALSE)





## get info on how many counties per state
births %>% 
  mutate(state = substr(county_name, nchar(county_name)-1, nchar(county_name))) %>%
  filter(year == "2016") %>%
  summarize(n = n(), .by = c("state"))

top_10_pct = plt_df %>% filter(total_cases_per_pop > 0, log_mean_birth_rate > log(min(plot_honeymoon$mu)), !is.na(mean_vax)) %>% 
  filter(total_cases_per_pop > quantile(total_cases_per_pop, 0.9)) 

top_10_pct %>% 
  mutate(min_tcpp = min(total_cases_per_pop), 
         max_tcpp = max(total_cases_per_pop)) %>%
  filter(total_cases_per_pop == max_tcpp)

mean_birth_rate = plt_df %>% filter(log_mean_birth_rate > log(min(plot_honeymoon$mu)), !is.na(mean_vax)) %>% mutate(birth_rate = exp(log_mean_birth_rate)) %>% 
  pull(birth_rate) %>% mean()

mean_vax_cov = plt_df %>% filter(log_mean_birth_rate > log(min(plot_honeymoon$mu)), !is.na(mean_vax)) %>%
  pull(mean_vax) %>% mean()

top_10_pct %>%
  mutate(above_mean_birth = exp(log_mean_birth_rate) > mean_birth_rate, 
         below_mean_vax = mean_vax < mean_vax_cov) %>%
  summarize(n = n(), pct = n()/nrow(top_10_pct), .by = c("above_mean_birth", "below_mean_vax"))
  
top_10_pct %>%
  mutate(above_mean_birth = exp(log_mean_birth_rate) > mean_birth_rate, 
         below_mean_vax = mean_vax < mean_vax_cov) %>% filter(!above_mean_birth, !below_mean_vax)

# counties with at least one case and vax cov > 0.95
plt_df %>% filter(total_cases_per_pop > 0, log_mean_birth_rate > log(min(plot_honeymoon$mu)), !is.na(mean_vax)) %>%
  mutate(tot = n()) %>%
  filter(mean_vax > 0.95, total_cases > 0) %>%
  summarize(n = n(), tot = mean(tot), pct = n()/mean(tot))

# counties with vax < 0.95
plt_df %>% filter(log_mean_birth_rate > log(min(plot_honeymoon$mu)), !is.na(mean_vax)) %>%
  mutate(vax_cat = ifelse(mean_vax < 0.90, "below 90%", ifelse(mean_vax > 0.95, "above 95%", "90%-95%"))) %>%
  summarize(n = n(), .by = "vax_cat") %>%
  mutate(pct = n/sum(n))

# counties with no cases and vax cov < 0.95
plt_df %>% filter(log_mean_birth_rate > log(min(plot_honeymoon$mu)), !is.na(mean_vax), mean_vax < 0.9) %>%
  mutate(tot = n()) %>%
  filter(total_cases != 0) %>%
  summarize(n = n(), tot = mean(tot), pct = n()/mean(tot))

