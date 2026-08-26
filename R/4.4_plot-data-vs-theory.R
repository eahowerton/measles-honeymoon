library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidytable)
library(cowplot)
library(dplyr)
library(beepr)
library(fields)
library(rootSolve)
library(patchwork)
library(stringr)
library(readxl)
library(metR)

### LOAD DATA ------------------------------------------------------------------
honeymoon_period_full = readRDS("data/output-data/honeymoon_period_by_birthrate.rds")
who_vacc = readRDS("data/output-data/WHO_vacc.rda")
who_drops_by_country = readRDS("data/output-data/drops_by_country_WHO.rda")
who_drops_by_country_summ = readRDS("data/output-data/drops_by_country_WHO_summary.rda")
us_vacc = readRDS("data/output-data/US_vacc.rda")
us_drops_by_county = readRDS("data/output-data/drops_by_county_US.rda")
us_drops_by_county_summ = readRDS("data/output-data/drops_by_county_US_summary.rda")
us_data_corrected = readRDS("data/output-data/us_data_corrected.rda")

rt_after_release_full_long = readRDS("data/output-data/rt_after_release_long.rds")
honeymoon_period_onemu = rt_after_release_full_long %>%
  filter(Rt > 1) %>%
  mutate(min_time = min(time), .by = c("waifw_id", "start_vax", "release_vax")) %>%
  filter(time == min_time) %>%
  select(-min_time)

### PLOT COUNTRY-LEVEL WHO VAX DROP RESULTS ------------------------------------
example_countries = c("Sudan", "Brazil", "Central African Republic", "Samoa", "Australia", "Benin")
start_yr = 2000
bin_width = 0.02
drop_bins = seq(-1, 0, bin_width)


# plot WHO drops by country
p_who_drops = ggplot(data = who_drops_by_country_summ %>% 
                       filter(nyears_data > 1) %>%
                       summarize(n = n(), .by = c("nyears_drop", "drop_bin")), 
            aes(x = nyears_drop, y = drop_bin - bin_width/2)) + # subtract bin_width/2 to get midpoint of bin on x-axis
  geom_line(data = honeymoon_period_onemu %>% filter(waifw_id == 5, start_vax == 0.95), 
            aes(x = time, y = -(1-release_vax)), linewidth = 0.4, linetype = "dotted") +
  geom_tile(aes(alpha = n), fill = "blue") +
  geom_point(data = who_drops_by_country %>% filter(country_name %in% example_countries),
             aes(y = drop), size = 1) +
  # label points, vary text positioning based on location of point to keep within plot window
  geom_text(data = who_drops_by_country %>%
              filter(country_name %in% setdiff(example_countries, c("Central African Republic", "Australia"))),
            aes(y = drop, label = paste0("\n", country_name)), size = 2) +
  geom_text(data = who_drops_by_country %>% filter(country_name == "Australia"),
            aes(y = drop, label = paste0("\n", country_name)), size = 2, hjust = 0) +
  geom_text(data = who_drops_by_country %>% filter(country_name == "Central African Republic"),
            aes(y = drop, label = "Central\nAfrican\nRepublic"),
            size = 2, hjust = 1, nudge_x = -0.3, vjust = 0.55) +
  scale_alpha_continuous(breaks = seq(2, 12, 2), name = "# countries") +
  scale_x_continuous(expand = c(0,0), name = "consecutive years dropping",
                     breaks = seq(0, 8, 2), limits = c(0.5, 8.4)) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent,
                     name = "largest coverage drop", limits = c(-1, 0)) +
  guides(alpha = guide_legend(title.position = "top", label.position = "bottom", nrow = 1,
                              keywidth = unit(0.3, "cm"), keyheight = unit(0.22, "cm"))) +
  theme_bw(base_size = 7) +
  theme(legend.position = c(0.68, 0.18),
        legend.direction = "horizontal",
        legend.key.spacing.x = unit(0, "pt"),
        panel.grid = element_blank())

# plot example coverage trajectories
p_who_exs = who_vacc %>% filter(country_name %in% example_countries, year > 1980) %>% 
  left_join(who_drops_by_country) %>% 
  mutate(drop_flag = ifelse(year >= start_yr & year <= start_yr + nyears_drop, TRUE, FALSE)) %>%
  filter(drop_flag, year > 2000) %>%
  ggplot(aes(x = year, y = coverage)) + 
  geom_line(data = who_vacc %>% filter(country_name %in% example_countries, year > start_yr), size = 0.4) +
  geom_line(color = "red", size = 0.6) +
  facet_wrap(vars(country_name), ncol = 2) + 
  scale_y_continuous(name = "MCV1 coverage", labels = scales::percent) + 
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank())


#### PLOT COUNTY-LEVEL US VAX DROP RESULTS -------------------------------------
example_FIPS = c("42003", "4013", "51680", "8031")
example_FIPS_labs = c("Maricopa, Arizona",
                      "Denver, Colorado",
                      "Allegheny, Pennsylvania",
                      "Lynchburg, Virginia" )
names(example_FIPS_labs) = c("4013", "8031", "42003", "51680")

p_us_drops = ggplot(data = us_drops_by_county_summ %>% summarize(n = n(), .by = c("nyears_drop", "drop_bin", "nyears_data")) %>%
              filter(nyears_data > 1),
            aes(x = nyears_drop, y = drop_bin-bin_width/2)) +
  geom_line(data = honeymoon_period_onemu %>% filter(waifw_id == 5, start_vax == 0.95), 
            aes(x = time, y = -(1-release_vax)), linewidth = 0.4, linetype = "dotted") +
  geom_tile(aes(alpha = n), fill = "purple") +
  geom_point(data = us_drops_by_county %>% filter(location_id %in% example_FIPS),
             aes(y = drop), size = 1) +
  # label points, vary text positioning based on location of point to keep within plot window
  geom_text(data = us_drops_by_county %>% filter(location_id %in% c("4013", "8031")),
            aes(y = drop, label = paste0("\n", county_name)), size = 2) +
  geom_text(data = us_drops_by_county %>% filter(location_id %in% c("42003", "51680")),
            aes(y = drop, label = county_name), size = 2, hjust = 0, nudge_x = 0.15) +
  scale_alpha_continuous(name = "# counties") +
  scale_x_continuous(expand = c(0,0), name = "consecutive years dropping", limits = c(0.5, 8)) +
  scale_y_continuous(expand = c(0,0), limits = c(-1, 0), labels = scales::percent,
                     name = "largest coverage drop") +
  guides(alpha = guide_legend(title.position = "top", label.position = "bottom", nrow = 1,
                              keywidth = unit(0.3, "cm"), keyheight = unit(0.22, "cm"))) +
  theme_bw(base_size = 7) +
  theme(legend.position = c(0.73, 0.18),
        legend.direction = "horizontal",
        legend.key.spacing.x = unit(0, "pt"),
        panel.grid = element_blank())

p_us_exs = us_vacc %>% filter(location_id %in% c(example_FIPS)) %>% 
  mutate(location_id = factor(location_id, levels = example_FIPS)) %>%
  left_join(us_drops_by_county_summ %>% mutate(location_id = as.factor(location_id))) %>%
  mutate(drop_flag = ifelse(start_year >= start_yr & start_year <= start_yr + nyears_drop, TRUE, FALSE)) %>%
  filter(drop_flag) %>%
  ggplot(aes(x = start_year,  y = value)) + 
  geom_line(data = us_vacc %>% filter(location_id %in% c(example_FIPS)) %>% 
              mutate(location_id = factor(location_id, levels = example_FIPS)), size = 0.4) + 
  geom_line(color = "red", size = 0.6) + 
  facet_wrap(vars(location_id), ncol = 2, labeller = labeller(location_id = example_FIPS_labs)) + 
  scale_x_continuous(name = "school year") +
  scale_y_continuous(name = "MMR coverage", limits = c(0, 1), labels = scales::percent) +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

#### PLOT US OUTBREKA DATA DATA ------------------------------------------------
tst_mu = seq(-5.5, -3.5, length.out = 20)
release_vax_full = seq(0.5, 1, 0.02)
plot_honeymoon = expand.grid(mu = exp(tst_mu), 
                             release_vax = release_vax_full) %>%
  left_join(honeymoon_period_full %>% select(mu, release_vax, time))

us_data_corrected %>%
  filter(births > 0) %>%
  pull(log_mean_birth_rate) %>%
  range(na.rm = TRUE)

# manual label placement for the highlighted counties 
# xb = births per 10,000 for the label anchor; y_lab = coverage for the label anchor.
# lab_color is black on the light top band and white over the dark region.
# point coordinates (px, py) are joined from the data so labels always connect to the true points.
e_label_pos = data.frame(
  location_name = c("Borden, Texas", "Spartanburg, South Carolina", "Terry, Texas",
                    "Yoakum, Texas", "Cochran, Texas", "Dawson, Texas", "Pawnee, Kansas",
                    "Mohave, Arizona", "Stevens, Kansas", "Gray, Kansas", "Kiowa, Kansas",
                    "Haskell, Kansas", "Gaines, Texas"),
  xb        = c(60, 56, 180, 215, 235, 150, 205, 52, 72, 180, 100, 140, 228),
  y_lab     = c(0.975, 0.935, 0.990, 0.955, 0.900, 0.892, 0.845, 0.858, 0.780,
                0.695, 0.600, 0.550, 0.730),
  hjust     = c(1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0),
  lab_color = c("black", "black", "black", "black", "white", "white", "white",
                "white", "white", "white", "white", "white", "white"),
  stringsAsFactors = FALSE
) %>%
  mutate(x_lab = log(xb / 1e4),
         label = gsub(", ", ",\n", location_name)) %>%
  left_join(us_data_corrected %>%
              transmute(location_name, px = log_mean_birth_rate, py = as.double(mean_vax)),
            by = "location_name")

p_us_outbreaks = ggplot(us_data_corrected %>%
                          filter(log_mean_birth_rate > log(min(plot_honeymoon$mu)))) + 
  geom_tile(data = plot_honeymoon, aes(x = log(mu), y = release_vax, fill = time)) +
  # geom_contour(data = plot_honeymoon, aes(x = log(mu), y = release_vax, z = time),
  #              color = "gray", breaks = c(3, 5, 7), linewidth = 0.2) +
  # metR::geom_text_contour(data = plot_honeymoon,
  #                         aes(x = log(mu), y = release_vax, z = time),
  #                         breaks = c(3, 5, 7), size = 2, color = "black",
  #                         stroke.colour = "gray",  # Outline color
  #                         stroke = 0.1) +
  geom_point(aes(x = log_mean_birth_rate, y = as.double(mean_vax)), 
             alpha = 0.8, shape = 21, color = "darkgray", size = 0.5, stroke = 0.2) +
  geom_segment(data = e_label_pos %>% filter(lab_color == "black"),
               aes(x = x_lab, y = y_lab, xend = px, yend = py),
               color = "black", linewidth = 0.3) +
  geom_segment(data = e_label_pos %>% filter(lab_color == "white"),
               aes(x = x_lab, y = y_lab, xend = px, yend = py),
               color = "white", linewidth = 0.3) +
  geom_text(data = e_label_pos %>% filter(lab_color == "black"),
            aes(x = x_lab, y = y_lab, label = label, hjust = hjust),
            color = "black", size = 2, lineheight = 0.85) +
  geom_text(data = e_label_pos %>% filter(lab_color == "white"),
            aes(x = x_lab, y = y_lab, label = label, hjust = hjust),
            color = "white", size = 2, lineheight = 0.85) +
  geom_point(data = us_data_corrected %>% 
               filter(total_cases_per_pop > 0, log_mean_birth_rate > log(min(plot_honeymoon$mu))),
             aes(x = log_mean_birth_rate, y = as.double(mean_vax), 
                 size = total_cases_per_pop, color = log(total_cases_per_pop*1e4))) +
  guides(size = "none") +
  scale_color_viridis_c(breaks = log(c(1, 10, 100)), #c(-2, 0, 2, 4, 6), 
                        labels = c(1, 10, 100),
                        #trans = scales::pseudo_log_trans(sigma = 0.001), 
                        na.value = "darkgray", name = "cases per\n10,000") +
  scale_fill_viridis_c(option = "rocket", na.value = "#FAEBDDFF", name = "time to\nRe > 1") +
  scale_size_continuous(range = c(0.5,3)) +
  scale_x_continuous(expand = c(0,0),
                     breaks = c(log(c(50, 100, 200)/1e4)), 
                     labels = c(50, 100, 200),
                     name = "births per 10,000") +
  scale_y_continuous(expand = c(0,0), name = "vaccination coverage", 
                     labels = scales::percent) +
  theme_bw(base_size = 7) + 
  theme(legend.key.width = unit(0.5, "cm"),
        legend.position = "bottom")

#### COMBINE INTO A SINGLE FIGURE ----------------------------------------------
plot_grid(
  plot_grid(
    plot_grid(p_who_exs, p_who_drops, nrow = 1, labels = c("A", "C"), 
              align = "h", axis = "tb", 
              label_size = 10, rel_widths = c(0.55, 0.45)),
    plot_grid(p_us_exs, p_us_drops, nrow = 1, labels = c("B", "D"), 
              align = "h", axis = "tb", label_size = 10, 
              rel_widths = c(0.55, 0.45)), 
    ncol = 1
  ),  
  p_us_outbreaks, labels = c(NA, "E"), rel_widths = c(0.55, 0.45), label_size = 10
  )

ggsave("figures/empirical_vax_declines.pdf", width = 8, height = 4)

