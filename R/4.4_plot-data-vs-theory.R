
honeymoon_period = saveRDS("data/output-data/honeymoon_period_by_birthrate.rds")
who_drops_by_country_summ = saveRDS("data/output-data/drops_by_country_WHO.rda")
drops_by_county_summ = saveRDS("data/output-data/drops_by_county_US.rda")

### PLOT COVERAGE DROPS --------------------------------------------------------

#### PLOT US DATA --------------------------------------------------------------
plot_honeymoon = expand.grid(mu = exp(tst_mu), 
                             release_vax = release_vax_full) %>%
  left_join(honeymoon_period %>% select(mu, release_vax, time))

plt_df %>% 
  filter(births > 0) %>%
  pull(log_mean_birth_rate) %>%
  range(na.rm = TRUE)

us_dat_fig = ggplot(plt_df %>% filter(log_mean_birth_rate > log(min(plot_honeymoon$mu)))) + 
  geom_tile(data = plot_honeymoon, aes(x = log(mu), y = release_vax, fill = time)) +
  geom_contour(data = plot_honeymoon, aes(x = log(mu), y = release_vax, z = time),
               color = "gray", breaks = c(3, 5, 7), linewidth = 0.2) +
  metR::geom_text_contour(data = plot_honeymoon,
                          aes(x = log(mu), y = release_vax, z = time),
                          breaks = c(3, 5, 7), size = 2, color = "black",
                          stroke.colour = "gray",  # Outline color
                          stroke = 0.1) +
  geom_point(aes(x = log_mean_birth_rate, y = as.double(mean_vax)), alpha = 0.8, shape = 21, color = "darkgray", size = 0.5, stroke = 0.2) +
  geom_point(data = plt_df %>% filter(total_cases_per_pop > 0, log_mean_birth_rate > log(min(plot_honeymoon$mu))),
             aes(x = log_mean_birth_rate, y = as.double(mean_vax), size = total_cases_per_pop, color = log(total_cases_per_pop*1e4))) +
  ggrepel::geom_text_repel(data = plt_df %>% filter(total_cases_per_pop > 0, log_mean_birth_rate > log(min(plot_honeymoon$mu))) %>% 
                             filter(total_cases_per_pop > quantile(total_cases_per_pop, 0.95)) %>% 
                             mutate(location_name2 = gsub(", ", ",\n",location_name)),
                           aes(x = log_mean_birth_rate, y = as.double(mean_vax), label = location_name2),
                           box.padding = unit(1.1, "lines"),
                           point.padding = unit(0.4, "lines"),
                           min.segment.length = 0,
                           max.overlaps = Inf,
                           segment.size = 0.4,
                           segment.color = 'black', size = 2
  ) +
  # geom_vline(xintercept = log(mean_birth_rate)) +
  # geom_hline(yintercept = mean_vax_cov) +
  guides(size = "none") +
  scale_color_viridis_c(breaks = log(c(1, 10, 100)), #c(-2, 0, 2, 4, 6), 
                        labels = c(1, 10, 100),
                        #trans = scales::pseudo_log_trans(sigma = 0.001), 
                        na.value = "darkgray", name = "cases per\n10,000") +
  scale_fill_viridis_c(option = "rocket", na.value = "#FAEBDDFF", name = "theoretical\nhoneymoon time") +
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
us_dat_fig

save("us_dat_fig",
     file="R/FINAL/data/us_dat_fig.rda")