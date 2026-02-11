library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(readxl)

folder = "data/US-data/"

# CA births from https://data.chhs.ca.gov/dataset/bc59b3db-3936-4c9b-9275-c5f4b3dc9023/resource/260fbfb0-e386-465d-85a4-6f28868dd51a/download/20241030_births_final_county_year_sup.csv 
# CA pop 2010-2020 https://dof.ca.gov/forecasting/demographics/estimates/annual-intercensal-population-estimates-by-race-ethnicity-with-age-and-gender-detail-2010-2020/
# US census county pop estimates 2020-2024 https://www.census.gov/data/tables/time-series/demo/popest/2020s-counties-total.html
# FIPS: https://data.transportation.gov/Railroads/State-County-and-City-FIPS-Reference-Table/eek5-pv8d/about_data

FIPS = read.csv(paste0(folder, "State__County_and_City_FIPS_Reference_Table_20251007.csv"))
FIPS_CA = FIPS %>% 
  filter(State.Name == "CALIFORNIA") %>%
  select(County.Name, County.Code) %>%
  unique()


vacc = read.csv(paste0(folder, "mmr_data_us_counties.csv")) %>% 
  melt(c("FIPS", "County", "State")) %>% 
  mutate(start_year = as.integer(substr(variable, 3,6)), 
         value = as.double(value)) #%>%
  # filter(State == "CA") # only CA for now where we have data
FIPS = vacc %>% 
  select(FIPS, County, State) %>% 
  unique()
births = read.csv(paste0(folder, "CA_births.csv")) %>%
  filter(Strata_Name == "Total Population", Geography_Type == "Residence") %>%
  left_join(FIPS %>% filter(State == "CA"))
pop1 = read_xlsx(paste0(folder, "CA_pop-2010-2020.xlsx"), sheet = 3, range = cell_cols("A:P")) %>%
  melt(c("FIPS", "Sex", "Race/ethnicity recode", "Age (0-100+)"), variable.name = "year") %>%
  summarize(value = sum(value), .by = c("FIPS", "year")) %>%
  mutate(year = as.integer(substr(year, nchar(as.character(year))-3, nchar(as.character(year))))) %>% 
  left_join(FIPS %>% filter(State == "CA"))
pop2 = read_xlsx(paste0(folder, "co-est2024-pop-06.xlsx"), skip = 4, n_max = 59, col_names = c("County", "base_estimate", "est_2020", "est_2021", "est_2022", "est_2023", "est_2024")) %>%
  select(-base_estimate) %>% 
  melt(c("County"), variable.name = "year") %>%
  filter(County != "California") %>%
  mutate(first_str = unlist(gregexpr(" ", as.character(County)))[1], .by = c("County")) %>%
  mutate(year = as.integer(substr(year, 5,8)), 
         County = substr(as.character(County), 2, first_str-1)) %>%
  select(-first_str) %>%
  left_join(FIPS %>% filter(State == "CA"))
pop = bind_rows(pop1 %>% mutate(source = "CA gov"), 
                pop2 %>% mutate(source = "US Census Bureau"))
pop = pop %>% filter(paste0(year, source) != "2020CA gov")

ggplot(data = vacc %>% filter(State == "CA") ) + 
  geom_line(aes(x = start_year, y = value, group = County), alpha = 0.2) + 
  theme_bw()

ggplot(data = pop %>% filter(!is.na(FIPS)), aes(x = year, y = value, color = source)) + 
  geom_line() + 
  geom_point() +
  facet_wrap(vars(FIPS), scales = "free")

ggplot(data = births %>% filter(Year > 2016), aes(x = Year, y = Count)) + 
  geom_line() + 
  geom_point() +
  facet_wrap(vars(FIPS), scales = "free")

combined = left_join(
  births %>% select(FIPS, Year, County, State, Count) %>% rename(births = Count, year = Year, county = County, state = State) %>% filter(year > 2016), 
  vacc %>% rename(county = County, state = State, year = start_year, coverage = value) %>% select(-variable)
  ) %>% 
  left_join(
    pop %>% select(-source) %>% rename(pop = value, county = County, state = State) %>% filter(year > 2016)
  )

# interpolate coverage where necessary
interp_vals = function(vals, yright = 1, yleft = 1){
  if(length(vals[!is.na(vals)]) < 3){return(rep(NA, length(vals)))}
  x = 1:length(vals)
  which_notna = which(!is.na(vals))
  return(approx(x[which_notna], vals[which_notna], x, yright = yright, yleft = yleft)$y)
}

combined = combined %>%
  mutate(interp_coverage = interp_vals(coverage), 
         interp_births = interp_vals(births, yleft = mean(births, na.rm = TRUE), yright = mean(births, na.rm = TRUE)),
         .by = c("county")) %>%
  # remove values > 1
  mutate(interp_coverage = ifelse(interp_coverage > 1, 1, interp_coverage)) %>%
  # reconstruct susceptibles
  mutate(new_susc = (1-interp_coverage)*births) 

combined_long = combined %>%
  melt(c("county", "year")) %>%
  mutate(value = as.double(value), year = as.double(year)) 

case_study_county = "Alameda"

combined_long %>% filter(county == case_study_county) %>%
  ggplot(aes(x = year, y = value)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(rows = vars(variable), scales = "free")

combined %>% 
  left_join(combined %>% 
              filter(!is.na(coverage)) %>% 
              summarize(n_coverage_years = n(), 
                        mean_cov = mean(coverage), .by = c("county"))) %>%
  ggplot(aes(x = year, y = coverage)) + 
  geom_hline(yintercept = 1) + 
  geom_line() + 
  geom_point() + 
  facet_wrap(vars(county))

# SOME SIMULATION OF SUSCEPTIBLE ACCUMULATION ----------------------------------
source("R/age-structure_deterministic/0_helper-functions.R")
# source("R/age-structure_deterministic/0_SIR-age.R")
# source("R/age-structure_deterministic/1_setup-WAIFW.R")

age_classes_data = 1:7
year_range = 2017:2024 # no birth rates in 2023

combined_filter = combined %>% 
  left_join(combined %>% 
              filter(!is.na(coverage)) %>% 
              summarize(n_coverage_years = n(), 
                        mean_cov = mean(coverage, na.rm = TRUE), .by = c("county")))

tmp_S = vector("list", length(year_range))
for(i in 1:length(year_range)){
  tmp_S[[i]] = combined_filter %>%
    filter(year <= year_range[i]) %>%
    rename(birth_year = year) %>%
    mutate(age = year_range[i] - birth_year + 1, 
           year = year_range[i])
}

tmp_S = bind_rows(tmp_S) %>%
  select(county, year, birth_year, age, new_susc, births) %>% 
  mutate(pop = sum(births), .by = c("county", "year")) %>%
  mutate(pct_s = new_susc/pop) %>% 
  filter(year >= (min(year_range) + max(age_classes_data)),
         age <= max(age_classes_data)) %>%
  right_join(expand.grid(year = (min(year_range) + max(age_classes_data)):max(year_range), 
                         age = age_classes_data)) %>%
  mutate(pct_s2 = ifelse(is.na(pct_s), 0, pct_s)) %>% 
  arrange(county, year)

#### USE THIS TO ESTIMATE RT (A BIT WEIRD, BECAUSE WHERE TO STOP) --------------
# NOTE: ALL OF THIS ASSUMES NO INFECTIONS (SO WORST CASE SCENARIO OF SORTS)

waifw_data = get_waifws(age_classes_data, background = 0.001, rescale = TRUE)

R0 = 14
paras = c(mu = 1/50, beta1 = 0,
          gamma = 365/14, delta = 0, p = 0) # delta = 1e-4
paras["beta0"]  = (paras["gamma"] + paras["mu"])*R0

# use PLOYMOD for now
Rt = tmp_S %>%
  summarize(Rt = get_Rt(waifw_data[[5]], pct_s2, beta = paras["beta0"], gamma = paras["gamma"], mu = paras["mu"], N = 1), 
            .by = c("year", "country_name")) %>% 
  mutate(max_Rt = max(Rt), min_Rt = min(Rt), 
         Rt_at_end = ifelse(Rt[year == max(year_range)] > 1, "Rt > 1", "Rt < 1"),
         .by = c("country_name")) %>%
  mutate(Rt_category = ifelse(min_Rt > 1, "Rt always > 1", ifelse(max_Rt < 1, "Rt always < 1", "Rt crosses 1")))

# now plot and see what happens
ggplot(data = Rt, aes(x = year, y = Rt, group = country_name, color = Rt_category)) + 
  geom_line(alpha = 0.7)+
  geom_hline(yintercept = 1, color = "black") + 
  # geom_label(data = Rt %>% filter(Rt_category == "Rt crosses 1" & Rt_at_end == "Rt > 1", year == max(year)),
  # aes(label = country_name)) +
  facet_wrap(vars(Rt_category)) +
  scale_color_brewer(palette = "Dark2") + 
  theme_bw() + 
  theme(legend.position = "none", 
        strip.background = element_blank())

ggplot(data = Rt, aes(x = year, y = Rt, color = Rt_category)) + 
  geom_line(aes(group = country_name), alpha = 0.5)+
  geom_hline(yintercept = 1, color = "black", linewidth = 0.8, linetype = "dashed") + 
  geom_line(data = Rt %>% filter(Rt_category == "Rt crosses 1" & Rt_at_end == "Rt > 1"), aes(group = country_name), color = "black")+
  # ggrepel::geom_text_repel(data = Rt %>% filter(Rt_category == "Rt crosses 1" & Rt_at_end == "Rt > 1", year == max(year)),
  #            aes(label = country_name), color = "black", size = 3) +
  geom_text(data = Rt %>% filter(Rt_category == "Rt crosses 1" & Rt_at_end == "Rt > 1", year == max(year)) %>% mutate(n = seq_len(n())),
            aes(x = 2022, y = 12.5 - n*0.5, label = country_name), color = "black", size = 2.25, hjust = 1) +
  geom_text(data = Rt %>% filter(Rt_category == "Rt always < 1", year == max(year)) %>% mutate(n = seq_len(n())),
            aes(x = 2022, y = 12.5 - n*0.5, label = country_name), size = 2.25, hjust = 1) +
  geom_text(data = Rt %>% select(country_name, Rt_category) %>% unique() %>% summarize(n = paste0("n = ", n()), .by = c("Rt_category")), 
            aes(x = 2000, y = 12, label = n), size = 5, hjust = 0) + 
  facet_wrap(vars(Rt_category)) +
  scale_color_brewer(palette = "Dark2") + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("figures/Rt_all_countries.pdf", width = 8, height = 3.75)

change_countries = Rt %>% filter(Rt_category == "Rt crosses 1" & Rt_at_end == "Rt > 1")

p1 = ggplot(data = combined_filter %>% filter(country_name %in% unique(change_countries$country_name), year > 2000)) + 
  geom_line(aes(x = year, y = interp_coverage), color = "darkgray") +
  geom_point(aes(x = year, y = coverage)) +
  facet_grid(cols = vars(country_name)) +
  labs(y = "MCV1 coverage") + 
  scale_y_continuous(labels = scales::percent)+
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
p2 = ggplot(data = tmp_S %>% filter(country_name %in% unique(change_countries$country_name))) + 
  geom_tile(aes(x = year, y = age, fill = pct_s2)) + 
  facet_grid(cols = vars(country_name)) +
  scale_fill_viridis_c(name = "% S") + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) + 
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
p3 = ggplot(data = change_countries, aes(x = year, y = Rt)) + 
  geom_line()+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed") + 
  facet_grid(cols = vars(country_name)) +
  scale_color_brewer(palette = "Dark2") + 
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
cowplot::plot_grid(p1, p2, p3, ncol = 1, align = "v", axis = "lr")
ggsave("figures/Rt_exemplar_countries.pdf", width = 8, height = 6)


# let's just look at vacc data (not try to reconstruct S)

pop = read_xlsx(paste0(folder, "co-est2024-pop_ALL.xlsx"), skip = 4, n_max = 3149, col_names = c("County", "base_estimate", "est_2020", "est_2021", "est_2022", "est_2023", "est_2024")) %>%
  select(-base_estimate) %>% 
  melt(c("County"), variable.name = "year") %>%
  filter(County != "United States") %>%
  mutate(first_str = unlist(gregexpr(" ", as.character(County)))[1],
         last_str = unlist(gregexpr(" ", as.character(County)))[length(unlist(gregexpr(" ", as.character(County))))], .by = c("County")) %>%
  mutate(year = as.integer(substr(year, 5,8)), 
         State = substr(as.character(County), last_str+1, nchar(as.character(County))),
         County = substr(as.character(County), 2, first_str-1)
  ) %>%
  select(-first_str, -last_str) %>%
  # left_join(usdata::county %>% select(name, state) %>% mutate(name = ) %>% rename(County = name, State = state))
  left_join(FIPS %>% mutate(row_id = seq_len(n())) %>% mutate(State = state.name[which(State == state.abb)], .by = row_id)) %>%
  filter(year == 2024)




ggplot(data = vacc, aes(x = start_year, y = value, group = FIPS)) + 
  geom_line(alpha = 0.2) + 
  facet_wrap(vars(State))

get_max_drop = function(vacc){
  n_yrs = length(vacc[!is.na(vacc)])
  if(n_yrs == 1){
    return(data.frame(nyears_data = n_yrs,
                      drop = NA, 
                      nyears_drop = NA))
  }
  vacc_diff = vector("list", n_yrs-1)
  drop = rep(NA, n_yrs-1)
  for(i in 1:(n_yrs-1)){
    vacc_diff[[i]] = lead(vacc, i) - vacc
    if(all(is.na(vacc_diff[[i]]))){drop[i] = NA;}
    else{drop[i] = min(vacc_diff[[i]], na.rm = TRUE)}
  }
  if(all(is.na(drop))){
    biggest_drop = NA
    which_biggest_drop = NA
  }
  else{
    biggest_drop = min(drop, na.rm = TRUE)
    which_biggest_drop = which.min(drop)
  }
  return(data.frame(nyears_data = n_yrs,
                    drop = biggest_drop, 
                    nyears_drop = which_biggest_drop))
}

drops_by_county = vacc %>% 
  arrange(State, FIPS, start_year) %>%
  dplyr::reframe(get_max_drop(value), .by = c("County", "FIPS", "State"))

drop_bins = seq(-0.5, 0, 0.01)
drops_by_county_summ  = drops_by_county %>%
  mutate(row_id = seq_len(n())) %>%
  mutate(drop_bin = ifelse(!is.na(drop), drop_bins[min(which(drop < drop_bins))], NA), .by = c("row_id"))


ggplot(data = drops_by_county_summ %>% summarize(n = n(), .by = c("nyears_drop", "drop_bin", "nyears_data")), 
       aes(x = nyears_drop, y = drop_bin)) + 
  geom_tile(aes(alpha = n), fill = "blue") + 
  facet_wrap(vars(nyears_data)) +
  theme_bw() + 
  theme(panel.grid = element_blank())


example_FIPS = c("42003", "4013", "51680", "8031")
example_FIPS_tmp = FIPS %>% 
  left_join(pop %>% select(County, value, FIPS) %>% rename(pop = value)) %>%
  filter(FIPS %in% example_FIPS) %>%
  mutate(lab = paste0(County, ", ", State, " (population: ", ifelse(pop/1e6 < 1, paste0(round(pop/1e3, 1), "K)"), paste0(round(pop/1e6, 1), "M)"))))
example_FIPS_labs = example_FIPS_tmp %>% pull(lab)
names(example_FIPS_labs) = example_FIPS_tmp %>% pull(FIPS)

# filter to four or more years of data (out of 7)
p2 = ggplot(data = drops_by_county_summ %>% summarize(n = n(), .by = c("nyears_drop", "drop_bin", "nyears_data")) %>% 
              filter(nyears_data > 3), 
            aes(x = nyears_drop, y = drop_bin)) + 
  geom_tile(aes(alpha = n), fill = "blue") + 
  geom_point(data = drops_by_county %>% filter(FIPS %in% example_FIPS), 
             aes(y = drop)) +
  geom_text(data = drops_by_county %>% filter(FIPS %in% example_FIPS), 
             aes(y = drop, label = paste0("\n",County, ", ", State)), size = 3) +
  scale_x_continuous(expand = c(0,0), name = "consecutive years dropping") + 
  scale_y_continuous(expand = c(0,0), limits = c(-0.8, 0), 
                     name = "largest coverage drop") + 
  theme_bw() + 
  theme(legend.position = "non", panel.grid = element_blank())

p1 = left_join(drops_by_county, pop %>% select(County, value, FIPS) %>% rename(pop = value)) %>%
  ggplot(aes(x = pop, y = drop, color = nyears_drop)) + 
  geom_point(alpha = 0.6) + 
  scale_x_log10(name = "log(population)") + 
  scale_y_continuous(name = "largest coverage drop") +
  scale_color_distiller(palette = "YlOrBr", direction = 1, name = "consecutive\nyears dropping") +
  theme_bw() + 
  theme(legend.position = c(0.55, 0.10), 
        legend.direction = "horizontal")


p3 = ggplot(data = vacc %>% filter(FIPS %in% c(example_FIPS)) %>% mutate(FIPS = factor(FIPS, levels = example_FIPS)), 
       aes(x = start_year,  y = value)) + 
  geom_line() + 
  facet_wrap(vars(FIPS), ncol = 1, labeller = labeller(FIPS = example_FIPS_labs)) + 
  scale_x_continuous(name = "school year") +
  scale_y_continuous(labels = scales::percent, name = "MMR coverage") +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

cowplot::plot_grid(p2, p3, p1, nrow = 1, rel_widths = c(0.3, 0.2, 0.3))
ggsave("R/FINAL/figures/US_vax.pdf", width = 10, height = 4)


# now try merging with pop size
FIPS = vacc %>% 
  select(FIPS, County, State) %>% 
  unique()






# now try fitting a linear regression through the data and get slope
regression = vacc %>% mutate(nyears_data = length(value[!is.na(value)]), .by = c("FIPS", "State")) %>%
  filter(nyears_data > 4) %>%
  select(FIPS, County, State) %>%
  unique() %>%
  mutate(slope = NA, lwr = NA, upr = NA)
for(i in 1:nrow(regression)){
  vacc_filt = vacc %>% filter(FIPS == regression[i, "FIPS"], State == regression[i, "State"])
  lm_tmp = lm(value ~ start_year, data = vacc_filt)
  regression[, "slope"] = lm_tmp$coefficients["start_year"]
  regression[, c("lwr", "upr")] = confint(lm_tmp, level = 0.95)[2, ]
}

regression %>%
  mutate(sign_both = sign(lwr) * sign(upr))

