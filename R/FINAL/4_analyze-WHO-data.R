library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)

folder = "data/WHO-data/"

births = read.csv(paste0(folder, "BirthRatePer1000downloadDec2024.csv")) %>%
  melt(c("Series.Name", "Series.Code", "Country.Name", "Country.Code"), variable.name = "year", value.name = "birth_rate") %>%
  mutate(year = as.integer(substr(year, 2, 5)))
vacc = read.csv(paste0(folder, "measlesVaccCoverFirstDose.csv")) %>%
  filter(COVERAGE_CATEGORY == "ADMIN") %>%
  rename(year = YEAR, coverage = COVERAGE) %>% 
  mutate(year = as.integer(year))
pop = read.csv(paste0(folder, "PopdownloadDec2024.csv")) %>%
  melt(c("X"), variable.name = "year", value.name = "pop") %>%
  mutate(year = as.integer(substr(year, 2, 5)))

combined = left_join(
  births %>% filter(Series.Code == "SP.DYN.CBRT.IN") %>% rename(country_name = Country.Name) %>% select(country_name, year, birth_rate) %>% mutate(birth_rate = as.double(birth_rate)), 
  vacc %>% rename(country_name = NAME) %>% select(country_name, year, coverage) %>% mutate(coverage = coverage/100), by = join_by(country_name, year)) %>% 
  left_join(
  pop %>% rename(country_name = X) %>% mutate(pop = as.double(pop)), by = join_by(country_name, year)
)

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
  mutate(interp_coverage = ifelse(interp_coverage > 1, 1, interp_coverage)) %>%
  # reconstruct susceptibles
  mutate(n_births = pop/1e3*interp_births, 
         new_susc = (1-interp_coverage)*n_births) 

combined_long = combined %>%
  melt(c("country_name", "year")) %>%
  mutate(value = as.double(value), year = as.double(year)) 

# case_study_country = "Denmark"
# case_study_country = "Paraguay"
case_study_country = "Portugal"
# case_study_country = "Malawi"

combined_long %>% filter(country_name == case_study_country, year > 1960) %>%
  ggplot(aes(x = year, y = value)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(rows = vars(variable), scales = "free")

combined %>% 
  filter(year > 1980) %>%
  left_join(combined %>% 
              filter(!is.na(coverage)) %>% 
              summarize(n_coverage_years = n(), 
                     mean_cov = mean(coverage), .by = c("country_name"))) %>%
  filter(n_coverage_years > 25) %>%
  ggplot(aes(x = year, y = coverage)) + 
  geom_hline(yintercept = 1) + 
  geom_line() + 
  geom_point() + 
  facet_wrap(vars(country_name))

# SOME SIMULATION OF SUSCEPTIBLE ACCUMULATION ----------------------------------
source("R/age-structure_deterministic/0_helper-functions.R")
# source("R/age-structure_deterministic/0_SIR-age.R")
# source("R/age-structure_deterministic/1_setup-WAIFW.R")

age_classes_data = 1:20
year_range = 1980:2022 # no birth rates in 2023

combined_filter = combined %>% 
  filter(year > 1980) %>%
  left_join(combined %>% 
              filter(!is.na(coverage)) %>% 
              summarize(n_coverage_years = n(), 
                     mean_cov = mean(coverage, na.rm = TRUE), .by = c("country_name"))) %>%
  filter(n_coverage_years > 30)

tmp_S = vector("list", length(year_range))
for(i in 1:length(year_range)){
  tmp_S[[i]] = combined_filter %>%
    filter(year <= year_range[i]) %>%
    rename(birth_year = year) %>%
    mutate(age = year_range[i] - birth_year + 1, 
           year = year_range[i])
}

tmp_S = bind_rows(tmp_S) %>%
  select(country_name, year, birth_year, age, new_susc, n_births) %>% 
  mutate(pop = sum(n_births), .by = c("country_name", "year")) %>%
  mutate(pct_s = new_susc/pop) %>% 
  filter(year >= (min(year_range) + max(age_classes_data)), 
         age <= max(age_classes_data)) %>% 
  right_join(expand.grid(year = (min(year_range) + max(age_classes_data)):max(year_range), 
                        age = age_classes_data), by = join_by(year, age)) %>%
  mutate(pct_s2 = ifelse(is.na(pct_s), 0, pct_s)) %>% 
  arrange(country_name, year)

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

