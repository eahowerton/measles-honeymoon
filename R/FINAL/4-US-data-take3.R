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

ggplot(data = honeymoon_period, aes(x = log(mu), y = release_vax, fill = time)) + 
  geom_tile()

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
  left_join(FIPS)

# do some manual updating of the data to associate each with a county
# get county-level detail on SC outbreak form dph.sc.gov (pulled on Feb 16, 2026)
# cherokee county is reported as <5 cases, so use 2 here (as floor(median))
# https://dph.sc.gov/diseases-conditions/infectious-diseases/measles-rubeola/measles-dashboard
upstate_sc = data.frame(county_name = c("Spartanburg", "Greenville", "Anderson", "Cherokee"), 
                        state_name = rep("South Carolina", 4), 
                        count = c(904, 35, 6, 2), 
                        location_id = c(45083, 45045, 45007, 45021)) %>%
  mutate(location_name = paste0(county_name, ", ", state_name))

# utah - can't disaggregate by county, so instead, we will aggregate birth rates and vax rates
# map from https://files.epi.utah.gov/Utah%20measles%20dashboard.html then converted using UT county map
utah_conversion = c(data_name = c(rep("Southwest", 5), 
                                  rep("Southeast", 3), 
                                  rep("Central", 6), 
                                  rep("TriCounty", 3), 
                                  rep("Bear River", 3)
                                  ), 
                    county_name = c("Beaver", "Iron", "Garfield", "Kane", "Washington", 
                                    "Grand", "Emery", "Carbon", 
                                    "Juab", "Millard", "Sanpete", "Sevier", "Piute", "Wayne", 
                                    "Daggett", "Uintah", "Duchesne", 
                                    "Box Elder", "Cache", "Rich"
                                    ))


# now add to measles cases
measles_cases = measles_cases %>% 
  filter(location_name != "Upstate, South Carolina") %>%
  bind_rows(upstate_sc)

# get measles vax
measles_vax = read.csv("data/US-data/mmr_data_us_counties.csv") %>%
  rename(location_id = FIPS, 
         state_name = State, 
         county_name = County) %>%
  mutate(SY2023_24 = as.double(SY2023_24), 
         location_name = paste0(county_name, ", ", state_name)) %>%
  mutate(mean_vax = mean(c(SY2017_18, SY2018_19, SY2019_20, SY2020_21, SY2021_22, SY2023_24), na.rm = TRUE), .by = c("location_name", "location_id", "county_name"))

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
         birth_rate = as.integer(birth_rate))


#### COMBINE AND PLOT ----------------------------------------------------------
plt_df = measles_cases %>%
  summarize(total_cases = sum(count), .by = c("location_id", "location_name")) %>%
  full_join(
    birth_pop %>% 
      mutate(location_id = as.integer(location_id)) %>%
      select(location_id, total_pop, births, birth_rate)
  ) %>%
  full_join(measles_vax %>% select(location_id, mean_vax)) %>%
  mutate(total_cases = ifelse(is.na(total_cases), 0, total_cases), 
         total_cases_per_pop = total_cases/total_pop) %>%
  mutate(log_mean_birth_rate = log(births/total_pop)) %>%
  filter(!is.na(location_id))

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

