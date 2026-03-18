library(ggplot2)
library(dplyr)
library(reshape2)
library(readxl)


folder = "data/WHO-data/"

start_yr = 2000

who_regions = read_xlsx(paste0(folder, "un-agencies-region-classification-for-country.xlsx"))

who_vacc = read.csv(paste0(folder, "measlesVaccCoverFirstDose.csv")) %>%
  filter(COVERAGE_CATEGORY == "WUENIC") %>%
  rename(year = YEAR, coverage = COVERAGE) %>% 
  mutate(year = as.integer(year)) %>%
  rename(country_name = NAME) %>% 
  select(country_name, year, coverage, CODE) %>% 
  mutate(coverage = coverage/100) %>% 
  filter(year > start_yr) %>%
  arrange(year) %>%
  mutate(coverage_diff = c(NA, diff(coverage)),
         average_coverage = mean(coverage, na.rm = TRUE), .by = c("country_name"))  %>%
  full_join(who_regions %>% 
              rename(country_name = `UNSD Country name`, 
                     region = `WHO Region name2`)
  ) %>% 
  filter(!is.na(region), !is.na(CODE))

ggplot(data = who_vacc, 
       aes(x = year, y = reorder(CODE, average_coverage))) + 
  geom_tile(aes(fill = coverage)) + 
  # geom_point(data = who_vacc %>% filter(coverage_diff < -0.05) %>% 
  #              mutate(coverage_drop_cat = factor(ifelse(coverage_diff < 0, ifelse(coverage_diff > -0.05, "no change - 5% decrease", ifelse(coverage_diff > -0.1, "5-10% decrease", ifelse( coverage_diff > -0.25, "10%-25% decrease", ">25% decrease"))), "increase")), levels = c("no change - 5% decrease", "5-10% decrease", "10%-25% decrease", ">25% decrease")), 
  #            aes(color = coverage_drop_cat), size = 0.8) +
  geom_point(data = who_vacc %>% filter(coverage_diff < -0.05), aes(color = ">5% decrease"), size = 0.6) +
  # geom_text(data = who_vacc %>% select(CODE, average_coverage) %>% unique(), 
  #            aes(x = 2023.5, label = paste0(round(average_coverage*100), "%")), size = 2, hjust = 1) +
  labs(fill = "MCV1\ncoverage", size = "annual change\nin coverage", color = "") + 
  facet_wrap(vars(region), scales = "free") +
  # scale_color_viridis_c() +
  # scale_color_manual(values = c("blue", "black", "darkblue")) +
  scale_color_manual(values = RColorBrewer::brewer.pal(7, "Greys")[rev(c(2,4,6)+1)]) +
  # scale_color_brewer(palette = "Greys", direction = -1) +
  scale_fill_distiller(palette = "OrRd", limits = c(0.15, 1)) +
  scale_shape_manual(values = c(19, NA, NA)) + 
  # scale_size_continuous(range = c(2,0), labels = scales::percent) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme_bw(base_size = 7) + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 4),
        panel.grid = element_blank(), 
        strip.background = element_blank())

ggsave("R/FINAL/figures/vacc_coverage_heatmap.pdf", width = 7.25, height = 5.25)