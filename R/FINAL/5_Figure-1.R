library(ggplot2)
library(dplyr)
library(reshape2)
library(readxl)
library(cowplot)

#### SETUP ---------------------------------------------------------------------
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
who_births = read.csv(paste0(folder, "BirthRatePer1000downloadDec2024.csv")) %>%
  melt(c("Series.Name", "Series.Code", "Country.Name", "Country.Code"), variable.name = "year", value.name = "birth_rate") %>%
  mutate(year = as.integer(substr(year, 2, 5)))
who_vacc = read.csv(paste0(folder, "measlesVaccCoverFirstDose.csv")) %>%
  filter(COVERAGE_CATEGORY == "WUENIC") %>%
  rename(year = YEAR, coverage = COVERAGE) %>% 
  mutate(year = as.integer(year))
who_pop = read.csv(paste0(folder, "PopdownloadDec2024.csv")) %>%
  melt(c("X"), variable.name = "year", value.name = "pop") %>%
  mutate(year = as.integer(substr(year, 2, 5)))
who_regions = read.csv(paste0(folder, "map.countries.regions.csv"))

combined = full_join(
  who_births %>% filter(Series.Code == "SP.DYN.CBRT.IN") %>% rename(country_name = Country.Name) %>% select(country_name, year, birth_rate) %>% mutate(birth_rate = as.double(birth_rate)), 
  who_vacc %>% rename(country_name = NAME) %>% select(country_name, year, coverage) %>% mutate(coverage = coverage/100)) %>% 
  full_join(
    who_pop %>% rename(country_name = X) %>% mutate(pop = as.double(pop))
  ) %>% 
  full_join(who_regions %>% rename(country_name = `Country...Region`) %>% melt(c("country_name"), variable.name = "region") %>% filter(!is.na(value)) %>% select(-value) %>% mutate(region = gsub(".region", "", region)))

# interpolate coverage where necessary
interp_vals = function(vals, yright = 1, yleft = 1){
  if(all(is.na(vals))){return(rep(NA, length(vals)))}
  x = 1:length(vals)
  which_notna = which(!is.na(vals))
  return(approx(x[which_notna], vals[which_notna], x, yright = yright, yleft = yleft)$y)
}

combined = combined %>%
  mutate(interp_coverage = interp_vals(coverage), 
         interp_births = interp_vals(birth_rate, yleft = mean(birth_rate, na.rm = TRUE), yright = mean(birth_rate, na.rm = TRUE)),
         .by = c("country_name")) %>%
  # remove values > 1
  mutate(interp_coverage = ifelse(interp_coverage > 1, 1, interp_coverage))

start_yr = 2000
who_drops_by_country = combined %>% 
  filter(year >= start_yr) %>%
  arrange(country_name, year) %>%
  dplyr::reframe(get_max_drop(coverage, start_yr), .by = c("country_name"))

bin_width = 0.02
drop_bins = seq(-1, 0, bin_width)
who_drops_by_country_summ  = who_drops_by_country %>%
  mutate(row_id = seq_len(n())) %>%
  mutate(drop_bin = ifelse(!is.na(drop), drop_bins[min(which(drop < drop_bins))], NA), .by = c("row_id"))

example_countries = c("Sudan", "Brazil", "Central African Republic", "Samoa", "Australia", "Benin") # "Burkina Faso", 

rt_after_release_full_long = readRDS("R/FINAL/data/rt_after_release_long.rds")
honeymoon_period = rt_after_release_full_long %>%
  filter(Rt > 1) %>%
  mutate(min_time = min(time), .by = c("waifw_id", "start_vax", "release_vax")) %>%
  filter(time == min_time) %>%
  select(-min_time)



# CHECK THESE REGION ASSIGNMENTS...
p0 = combined %>%
  filter(year > 1980, !is.na(region)) %>%
  reframe(quantile = paste0("Q", c(25, 50, 75)), 
            value = quantile(coverage, c(0.25, 0.5, 0.75), na.rm = TRUE), .by = c("year", "region")) %>%
  dcast(region + year ~ quantile, value.var = "value") %>% 
  ggplot(aes(x = year)) + 
  geom_ribbon(aes(ymin = Q25, ymax = Q75, fill = region), alpha = 0.2) + 
  geom_line(aes(y = Q50, color = region)) + 
  labs(y = "MCV1 coverage") +
  guides(color = guide_legend(ncol = 2), fill = guide_legend(ncol = 2)) + 
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Accent")[5:8],
                     labels = c("AFRO", "AMRO", "WPRO", "SEARO")) + 
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Accent")[5:8],
                    labels = c("AFRO", "AMRO", "WPRO", "SEARO")) + 
  theme_bw(base_size = 7) + 
  theme(legend.position = c(0.65, 0.2), 
        legend.direction = "horizontal",
        panel.grid = element_blank())
p1 = ggplot(data = who_drops_by_country_summ %>% filter(nyears_data > 1) %>% summarize(n = n(), .by = c("nyears_drop", "drop_bin")), 
       aes(x = nyears_drop, y = drop_bin - bin_width/2)) + # subtract bin_width/2 to get midpoint of bin on x-axis
  geom_line(data = honeymoon_period %>% filter(waifw_id == 5, start_vax == 0.95), 
            aes(x = time, y = -(1-release_vax)), linewidth = 0.4, linetype = "dotted") +
  geom_tile(aes(alpha = n), fill = "blue") +
  geom_point(data = who_drops_by_country %>% filter(country_name %in% example_countries),
             aes(y = drop), size = 1) +
  geom_text(data = who_drops_by_country %>% filter(country_name %in% example_countries),
            aes(y = drop, label = paste0("\n", country_name)), size = 2) +
  scale_x_continuous(expand = c(0,0), name = "consecutive years dropping",
                      breaks = seq(0, 8, 2)) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent,
                     name = "largest coverage drop", limits = c(-1, 0)) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none", 
        panel.grid = element_blank())
p2 = combined %>% filter(country_name %in% example_countries, year > 1980) %>% 
  left_join(who_drops_by_country) %>% 
  mutate(drop_flag = ifelse(year >= start_yr & year <= start_yr + nyears_drop, TRUE, FALSE)) %>%
  filter(drop_flag) %>%
  ggplot(aes(x = year, y = coverage)) + 
  geom_vline(xintercept = start_yr, linetype = "dotted", size = 0.4) + 
  geom_line(data = combined %>% filter(country_name %in% example_countries, year > 1980), size = 0.4) +
  geom_line(color = "red", size = 0.6) +
  facet_wrap(vars(country_name), ncol = 2) + 
  scale_y_continuous(name = "MCV1 coverage", labels = scales::percent) + 
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank())
p_who = plot_grid(p0, p1, p2, nrow = 1, labels = c("A", "B", "C"))

# finally plot the increases globally


 #### COVERAGE DROPS IS US DATA -------------------------------------------------
# US census county pop estimates 2020-2024 https://www.census.gov/data/tables/time-series/demo/popest/2020s-counties-total.html
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

pop = read_xlsx("data/US-data/co-est2024-pop_ALL.xlsx", skip = 4, n_max = 3149, col_names = c("County", "base_estimate", "est_2020", "est_2021", "est_2022", "est_2023", "est_2024")) %>%
  select(-base_estimate) %>% 
  melt(c("County"), variable.name = "year") %>%
  filter(County != "United States") %>%
  mutate(first_str = unlist(gregexpr(",", as.character(County)))[1], .by = c("County")) %>%
  mutate(year = as.integer(substr(year, 5,8)), 
         state_name = substr(as.character(County), first_str+2, nchar(as.character(County))),
         county_name = substr(as.character(County), 2, first_str-1)
  ) %>%
  mutate(county_name = gsub(" County", "", county_name)) %>%
  mutate(county_name = gsub(" Parish", "", county_name)) %>%
  mutate(county_name = gsub(" city", "", county_name)) %>%
  filter(year == 2024) %>%
  # select(-first_str, -County) %>%
  left_join(FIPS) 


drops_by_county = vacc %>% 
  arrange(state_abbrev, location_id, start_year) %>%
  dplyr::reframe(get_max_drop(value, 2017, 2023), .by = c("county_name", "location_id", "state_abbrev"))

drop_bins = seq(-0.5, 0, 0.01)
drops_by_county_summ  = drops_by_county %>%
  mutate(row_id = seq_len(n())) %>%
  mutate(drop_bin = ifelse(!is.na(drop), drop_bins[min(which(drop < drop_bins))], NA), .by = c("row_id"))

example_FIPS = c("42003", "4013", "51680", "8031")
example_FIPS_tmp = FIPS %>%
  filter(location_id %in% example_FIPS) %>%
  mutate(lab = paste0(county_name, ", ", state_name))
example_FIPS_labs = example_FIPS_tmp %>% pull(lab)
names(example_FIPS_labs) = example_FIPS_tmp %>% pull(location_id)
  
# example_FIPS_tmp = pop %>% rename(pop = value) %>%
#   filter(location_id %in% example_FIPS) %>%
#   mutate(lab = paste0(county_name, ", ", state_name, "\n(population: ", ifelse(pop/1e6 < 1, paste0(round(pop/1e3, 1), "K)"), paste0(round(pop/1e6, 1), "M)"))))
# example_FIPS_labs = example_FIPS_tmp %>% pull(lab)
# names(example_FIPS_labs) = example_FIPS_tmp %>% pull(location_id)

# filter to four or more years of data (out of 7)
p4 = left_join(drops_by_county, pop %>% select(county_name, value, location_id) %>% rename(pop = value) %>% filter(pop > 1e3)) %>%
  ggplot(aes(x = pop, y = drop, color = nyears_drop)) +
  geom_point(alpha = 0.6) +
  scale_x_log10(name = "population (log scale))", ) +
  scale_y_continuous(expand = c(0,0), limits = c(-1,0),
                     name = "largest coverage drop", labels = scales::percent) +
  scale_color_distiller(palette = "YlOrBr", direction = 1, name = "consecutive\nyears dropping") +
  theme_bw(base_size = 7) +
  theme(legend.position = c(0.6, 0.15),
        legend.direction = "horizontal", 
        panel.grid = element_blank())

p5 = ggplot(data = drops_by_county_summ %>% summarize(n = n(), .by = c("nyears_drop", "drop_bin", "nyears_data")) %>%
              filter(nyears_data > 1),
            aes(x = nyears_drop, y = drop_bin)) +
  geom_line(data = honeymoon_period %>% filter(waifw_id == 5, start_vax == 0.95), 
            aes(x = time, y = -(1-release_vax)), linewidth = 0.4, linetype = "dotted") +
  geom_tile(aes(alpha = n), fill = "purple") +
  geom_point(data = drops_by_county %>% filter(location_id %in% example_FIPS),
             aes(y = drop), size = 1) +
  geom_text(data = drops_by_county %>% filter(location_id %in% example_FIPS),
            aes(y = drop, label = paste0("\n", county_name, ", ", state_abbrev)), size = 2) +
  scale_x_continuous(expand = c(0,0), name = "consecutive years dropping", limits = c(0.5, 8)) +
  scale_y_continuous(expand = c(0,0), limits = c(-1, 0), labels = scales::percent,
                     name = "largest coverage drop") +
  theme_bw(base_size = 7) +
  theme(legend.position = "non", panel.grid = element_blank())

p6 = vacc %>% filter(location_id %in% c(example_FIPS)) %>% 
  mutate(location_id = factor(location_id, levels = example_FIPS)) %>%
  left_join(drops_by_county_summ %>% mutate(location_id = as.factor(location_id))) %>%
  mutate(drop_flag = ifelse(start_year >= start_yr & start_year <= start_yr + nyears_drop, TRUE, FALSE)) %>%
  filter(drop_flag) %>%
  ggplot(aes(x = start_year,  y = value)) + 
  geom_line(data = vacc %>% filter(location_id %in% c(example_FIPS)) %>% 
              mutate(location_id = factor(location_id, levels = example_FIPS)), size = 0.4) + 
  geom_line(color = "red", size = 0.6) + 
  facet_wrap(vars(location_id), ncol = 2, labeller = labeller(location_id = example_FIPS_labs)) + 
  scale_x_continuous(name = "school year") +
  scale_y_continuous(name = "MMR coverage", limits = c(0, 1), labels = scales::percent) +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

p_us = cowplot::plot_grid( p4, p5, p6, nrow = 1, labels = c("D", "E", "F"))

plot_grid(p_who, p_us, ncol = 1)

ggsave("R/FINAL/figures/empirical_vax_declines.pdf", width = 12, height = 6)

ggsave("R/FINAL/figures/US_vax.pdf", width = 10, height = 4)

plot_grid(p0, 
  plot_grid(
    plot_grid(p1, p2, nrow = 1, labels = c("B", "C")), 
    cowplot::plot_grid(p5, p6, nrow = 1, labels = c("D", "E")), ncol = 1
  ), labels = c("A", NA), nrow = 1, rel_widths = c(0.4, 0.6))

ggsave("R/FINAL/figures/empirical_vax_declines_v2.pdf", width = 12, height = 6)

## look for some outbreaks
who_cases <- read.csv(paste0(folder, "measles-cases-by-month.csv"),stringsAsFactors=FALSE) %>%
  rename(country_name = Country, 
         year = Year) %>% 
  melt(c("country_name", "Region", "ISO3", "year")) %>% 
  mutate(continuous_time = year + match(variable, month.name)/12) %>% 
  left_join(who_pop %>% rename(country_name = X) %>% mutate(pop = as.double(pop)) %>% 
              bind_rows(who_pop %>% filter(year == 2023) %>% rename(country_name = X) %>% mutate(pop = as.double(pop), year = 2024)))

example_countries2 = c("Brazi", "Zimbabwe", "Tunisia", "Vietnam", "Ukraine", "Krygyzstan", "Kazakstan", "Gabon", "Greece", "Cambodia", "Israel", "Serbia",  "Madagascar", "Equatorial Guinea") #Egypt","Slovakia",

q0 = ggplot(data = who_cases %>% filter(country_name %in% example_countries2)) + 
  geom_line(aes(x = continuous_time, y = value/pop * 1e4)) +
  facet_wrap(vars(country_name), ncol = 2, scales = "free") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank())
q1 = ggplot(data = who_drops_by_country_summ %>% filter(nyears_data > 1) %>% summarize(n = n(), .by = c("nyears_drop", "drop_bin")), 
            aes(x = nyears_drop, y = drop_bin - bin_width/2)) + # subtract bin_width/2 to get midpoint of bin on x-axis
  geom_tile(aes(alpha = n), fill = "blue") +
  geom_point(data = who_drops_by_country %>% filter(country_name %in% example_countries2),
             aes(y = drop)) +
  geom_text(data = who_drops_by_country %>% filter(country_name %in% example_countries2),
            aes(y = drop, label = paste0("\n", country_name)), size = 3) +
  scale_x_continuous(expand = c(0,0), name = "consecutive years dropping",
                     breaks = seq(0, 8, 2)) +
  scale_y_continuous(expand = c(0,0),
                     name = "largest coverage drop", limits = c(-1, 0)) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank())
q2 = combined %>% filter(country_name %in% example_countries2, year > 1980) %>% 
  left_join(who_drops_by_country) %>% 
  mutate(drop_flag = ifelse(year >= start_yr & year <= start_yr + nyears_drop, TRUE, FALSE)) %>%
  filter(drop_flag) %>%
  ggplot(aes(x = year, y = coverage)) + 
  geom_vline(xintercept = start_yr, linetype = "dotted") + 
  geom_line(data = combined %>% filter(country_name %in% example_countries2, year > 1980)) +
  geom_line(color = "red", size = 1) +
  facet_wrap(vars(country_name), ncol = 2) + 
  scale_y_continuous(name = "MCV1 coverage") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank())
plot_grid(q0, q1, q2, nrow = 1, labels = c("A", "B", "C"))
ggsave("R/FINAL/figures/outbreaks_and_declines.pdf", width = 12, height = 6)



### now some extra figures
ggplot(data = who_cases %>% filter(country_name %in% example_countries2)) + 
  geom_line(aes(x = continuous_time, y = value)) +
  geom_point(aes(x = continuous_time, y = value), size = 0.8) +
  facet_wrap(vars(country_name), scales = "free", ncol = 5) + 
  labs(y = "reported cases") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        panel.grid = element_blank(), 
        strip.background = element_blank())
ggsave("R/FINAL/figures/honeymoon_examples.pdf", width = 10, height = 4)


### ANOTHER VERSIO
load("R/FINAL/data/us_dat_fig.rda")

plot_grid(
  plot_grid(
  plot_grid(p2, p1, nrow = 1, labels = c("A", "C"), align = "h", axis = "tb", label_size = 10, rel_widths = c(0.55, 0.45)),
  plot_grid(p6, p5, nrow = 1, labels = c("B", "D"), align = "h", axis = "tb", label_size = 10, rel_widths = c(0.55, 0.45)), 
  ncol = 1
),  us_dat_fig, labels = c(NA, "E"), rel_widths = c(0.55, 0.45), label_size = 10)

ggsave("R/FINAL/figures/empirical_vax_declines_v3.pdf", width = 8, height = 4)


