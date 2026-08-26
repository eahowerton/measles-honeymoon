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
  geom_point(data = who_vacc %>% filter(coverage_diff < -0.05), aes(color = ">5% decrease"), size = 0.6) +
  labs(fill = "MCV1\ncoverage", size = "annual change\nin coverage", color = "") + 
  facet_wrap(vars(region), scales = "free") +
  scale_color_manual(values = RColorBrewer::brewer.pal(7, "Greys")[rev(c(2,4,6)+1)]) +
  scale_fill_distiller(palette = "OrRd", limits = c(0.15, 1)) +
  scale_shape_manual(values = c(19, NA, NA)) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) +
  theme_bw(base_size = 7) + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 4),
        panel.grid = element_blank(), 
        strip.background = element_blank())

ggsave("figures/vacc_coverage_heatmap.pdf", width = 7.25, height = 5.25)

