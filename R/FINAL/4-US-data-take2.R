# births_std = births %>%
#   rename(location_name_full = county_name) %>%
#   mutate(state_abbrev = substr(location_name_full, nchar(location_name_full)-1, nchar(location_name_full))) %>%
#   mutate(row_id = seq_len(n())) %>%
#   mutate(county_name = substr(location_name_full, 1, gregexpr(" County,", location_name_full)[[1]][1] - 1), 
#          state_name =  ifelse(state_abbrev == "DC", "District of Columbia", state.name[which(state.abb == state_abbrev)]), .by = c("row_id")) %>%
#   select(location_id, state_abbrev, state_name, county_name, births, total_pop, birth_rate, year) %>%
#   mutate(location_id = as.integer(location_id))

measles_cases_std = read.csv("data/US-data/measles_county_all_updates.csv") %>%
  mutate(county_name = substr(location_name, 1, gregexpr(",", location_name)[[1]][1] - 1), 
         state_name = substr(location_name, gregexpr(", ", location_name)[[1]][1] + 2, nchar(location_name)),
         # state_abbrev = ifelse(state_name == "District of Columbia", "DC", state.abb[which(state.name == state_name)]),
         .by = c("location_name", "location_id")) %>%
  left_join(data.frame(state_abbrev = c(state.abb, "DC"), 
                       state_name = c(state.name, "District of Columbia")))

summ_measles_cases_std = measles_cases_std %>%
  summarize(total_cases = sum(value), .by = c(state_abbrev, state_name, county_name))

birth_pop = read.csv("data/US-data/social_explorer_2019_US_census_birth_pop.csv", stringsAsFactors = FALSE)[-1,] %>%
  mutate(county_name = gsub(" County", "", Name.of.Area), 
         state_name = substr(Qualifying.Name, gregexpr(", ", Qualifying.Name)[[1]][1] + 2, nchar(Qualifying.Name)), 
         .by = c("Qualifying.Name")) %>% 
  rename(total_pop = `Total.Population`, 
         births = `Births`, 
         birth_rate = `Births.Rate.per.1.000.Population`, 
         location_id = FIPS) %>%
  select(county_name, state_name, location_id, total_pop, births, birth_rate)

vax_std = read.csv("data/US-data/mmr_data_us_counties.csv") %>% 
  melt(c("FIPS", "County", "State")) %>% 
  mutate(start_year = as.integer(substr(variable, 3,6)), 
         value = as.double(value)) %>%
  rename(county_name = County, 
         state_abbrev = State, 
         location_id = FIPS) %>% 
  left_join(data.frame(state_abbrev = c(state.abb, "DC"), 
                       state_name = c(state.name, "District of Columbia"))) %>%
  mutate(year = as.integer(substr(variable, 3, 6))) %>% 
  select(-variable) %>% 
  rename(vacc_pct = value) %>% 
  mutate(vacc_pct = as.double(vacc_pct))

combined_std = left_join(summ_measles_cases_std, 
                         birth_pop) %>% 
  mutate(location_id = as.integer(location_id)) %>%
  left_join(vax_std %>% mutate(location_id = as.integer(location_id)))

combined_avg = left_join(summ_measles_cases_std, 
                         birth_pop) %>% 
  mutate(location_id = as.integer(location_id)) %>%
  left_join(vax_std %>% summarize(vacc_pct = mean(vacc_pct, na.rm = TRUE), .by = c("location_id", "state_abbrev", "state_name", "county_name")) %>% mutate(location_id = as.integer(location_id))) %>%
  mutate(total_cases = as.integer(total_cases), 
         total_pop = as.integer(total_pop), 
         birth_rate = as.double(birth_rate))


ggplot(data = combined_avg, aes(x = birth_rate, y = vacc_pct)) + 
  geom_point(aes(size = total_cases/total_pop, color = total_cases/total_pop)) + 
  scale_color_viridis_c()

# now add a few from social explorer (2019)
social_explorer_manual = data.frame(county_name = c("Gaines", "Lea", "Terry", "Oconto", "Gray", "Lamar", "Dawson", 
                                                    "Haskell", "Yoakum", "Gallatin", "Williams", "Cochran", "Knox", "Luna", 
                                                    "Grand Forks", "Stevens", "Lincoln", "Carbon", "Dallam", "Hockley", "Osceola", 
                                                    "Pawnee"), 
           state_name = c("Texas", "New Mexico", "Texas", "Wisconsin", "Kansas", "Texas", "Texas", 
                          "Kansas", "Texas", "Montana", "North Dakota", "Texas", "Ohio", "New Mexico", 
                          "North Dakota", "Kansas", "South Dakota", "Wyoming", "Texas", "Texas", "Michigan", 
                          "Kansas"), 
           state_abbrev = c("TX", "NM", "TX", "WI", "KS", "TX", "TX",
                            "KS", "TX", "MT", "ND", "TX", "OH", "NM", 
                            "ND", "KS", "SD", "WY", "TX", "TX", "MI", 
                            "KS"), 
           births_man = c(430, 1086, 181, 369, 94, 619, 181, 
                      49, 151, 1262, 725, 41, 741, 346, 
                      974, 62, 866, 176, 151, 306, 265, 
                      76), 
           birth_rate_man = c(20.33, 15.45, 14.67, 9.74, 15.58, 12.43, 14.29, 
                          12.31, 17.45, 11.16, 19.85, 14.43, 11.93, 14.54, 
                          13.93, 11.25, 14.44, 11.86, 20.70, 13.31, 11.32, 
                          11.70), 
           total_pop_man = c(21492, 71070, 12337, 37930, 5988, 49859, 12728, 
                         3968, 8713, 114434, 37589, 2853, 62322, 23709, 
                         69451, 5485, 61128, 14800, 7287, 23021, 23460, 
                         6414)
           )


combined_se_supp = combined_avg #%>%  
  # left_join(social_explorer_manual) %>% 
  # mutate(births = ifelse(is.na(births), births_man, births), 
         # total_pop = ifelse(is.na(total_pop), total_pop_man, total_pop), 
         # birth_rate = ifelse(is.na(birth_rate), birth_rate_man, birth_rate)
         # )

ggplot(data = combined_se_supp, aes(x = birth_rate, y = vacc_pct)) + 
  geom_point(aes(size = log(total_cases), color = log(total_cases))) + 
  scale_color_viridis_c() +
  theme_bw()

## a nice one
broader_return_time = expand.grid(log_mu = seq(log(1/50), log(1/130), length.out = 12),
                                  # mu =seq(1/50, 1/130, length.out = 12), 
                                  beta0 = paras["beta0"], 
                                  v = 1, 
                                  gamma = paras["gamma"], 
                                  new_v = seq(0, 0.99, 0.01)) %>%
  mutate(mu = exp(log_mu), 
         R0 = (beta0)/((gamma + mu)), 
         herd_imm = 1-1/R0,
         S_star = ifelse(v > herd_imm, 1-v,
                         (mu + gamma)/beta0)) %>%
  mutate(time_to_cross = -1/mu * log(((1/R0) - (1-new_v))/(S_star - (1-new_v))))

ggplot(data = combined_se_supp, aes(x = birth_rate/1000, y = vacc_pct)) + 
  geom_tile(data = broader_return_time, aes(x = mu, y = new_v, fill = R0)) + 
  geom_contour(data = broader_return_time, aes(x = mu, y = new_v, z = R0), 
               color = "black", breaks = c(1, 3, 5, 10)) + 
  metR::geom_label_contour(data = broader_return_time, 
                           aes(x = mu, y = new_v, z = R0), 
                           breaks = c(1, 3, 5, 10)) + 
  geom_point(aes(size = total_cases/total_pop*1000, color = total_cases/total_pop*1000)) + 
  scale_color_viridis_c(trans = "log") +
  scale_fill_distiller(palette = "Greys") + 
  scale_size_continuous(trans = "log") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()

# another option where vacc pct is not required
ggplot(data = combined_se_supp, aes(x = birth_rate, y = total_cases/total_pop)) + 
  geom_point(aes(color = vacc_pct), size = 2.5) + 
  scale_color_viridis_c(trans = "log", labels = scales::percent, na.value = "gray", direction = -1, limits = c(0.75, 1)) + 
  scale_y_log10() + 
  theme_bw()

# another nice one
ggplot(data = combined_se_supp, aes(x = vacc_pct, y = total_cases/total_pop)) + 
  geom_point(aes(color = birth_rate), size = 2.5) +
  coord_cartesian(xlim = c(0.75, 1)) +
  scale_color_viridis_c(trans = "log") +
  scale_y_log10() +
  theme_bw()

ggplot(data = combined_se_supp %>% filter(vacc_pct > 0.6), aes(x = total_pop, y = birth_rate)) + 
  geom_point(aes(size = log(total_cases), color = vacc_pct)) + 
  scale_color_viridis_c() + 
  scale_x_log10()


# BACKGROUND?
broader_return_time = expand.grid(mu = seq(1/50, 1/80, length.out = 6), 
                                  beta0 = paras["beta0"], 
                                  v = pre_v, 
                                  gamma = paras["gamma"], 
                                  new_v = seq(0, 0.90, 0.01)) %>%
  mutate(R0 = (beta0)/((gamma + mu)), 
         herd_imm = 1-1/R0,
         S_star = ifelse(v > herd_imm, 1-v,
                         (mu + gamma)/beta0)) %>%
  mutate(time_to_cross = -1/mu * log(((1/R0) - (1-new_v))/(S_star - (1-new_v))))

ggplot(data = broader_return_time, aes(x = mu, y = new_v, fill = time_to_cross)) +
  geom_tile() + 
  geom_contour(aes(z = time_to_cross), color = "black", breaks = c(1, 3, 5, 10)) + 
  metR::geom_label_contour(aes(z = time_to_cross), breaks = c(1, 3, 5, 10)) + 
  facet_wrap(vars(v)) + 
  scale_fill_distiller(palette = "Greys") + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0))

