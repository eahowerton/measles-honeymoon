library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidytable)
library(epimdr2)
library(cowplot)
library(dplyr)
library(beepr)

source("R/age-struc-functions.R")

#### MODEL SETUP ---------------------------------------------------------------
compartments = c("S", "E", "I", "R")
times = seq(0, 100, by = 1/365)
times_long = seq(0, 300, by = 1/365)

no_vax = approxfun(times_long, rep(0, length(times_long)))

paras = c(mu = 1/50, N = 1, beta0 = 1000, beta1 = 0.2, omega = 0,
          sigma = 365/8, gamma = 365/5, vax_pct = 0, delta = 0, p = 0)

start_pop = 1
fert = paras["mu"]
mort = paras["mu"]

#### TEST WITH ONE COMPARTMENT -------------------------------------------------
age_classes1 =  c(100)
rslts1 <- run_ode(
  age_classes = age_classes1, mort = mort, fert = fert, start_pop = start_pop, 
  IC_type = "std", max_t = max(times_long), params = paras,
  plot_flag = TRUE, plot_title = "1 age class")

#### TEST WITH TWO COMPARTMENTS ------------------------------------------------
age_classes2 = c(20, 40)
rslts2 <- run_ode(
  age_classes = age_classes2, mort = rep(mort, length(age_classes2)),
  fert = rep(fert, length(age_classes2)), start_pop = start_pop,
  IC_type = "std", max_t = max(times_long)+100, params = paras,
  plot_flag = TRUE, plot_title = "2 age classes")

#### TEST WITH REALISTIC STRUCTURE ---------------------------------------------
# up to 5 years in months, up to 10 in years, and up to 70 in 5 years 
# (multiply by 4 to go from months to weeks)
age_classes = c(seq(0.5, 5, 1/2), 6:9, seq(10,70,by=10))
bin_width = diff(c(0, age_classes))
names(bin_width) = age_classes

start.time <- Sys.time()
rslts3 <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop,
  IC_type = "std", max_t = max(times_long), params = paras, beep_flag = TRUE,
  plot_flag = TRUE, plot_title = "realistic structure, constant WAIFW")
Sys.time() - start.time

#### NOW ADD VACCINATION -------------------------------------------------------
vax_change_times = c(0, 290, 295)
vax_rates = c(0, 0.8, 0.4)

start.time <- Sys.time()
rslts3_v <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop, 
  vax_change_times = vax_change_times, vax_rates = vax_rates,
  IC_type = "std", max_t = max(times_long), params = paras, beep_flag = TRUE,
  plot_flag = TRUE, plot_title = "realistic structure, constant WAIFW")
Sys.time() - start.time

#### TEST WITH POLYMOD CONTACTS ------------------------------------------------
# now add more realistic mixing with POLYMOD data
W = create_polymod_matrix(age_classes, plot_flag = TRUE, age_classes_to_label = c(1, 3, 5, 7, 9, 10, 30, 50, 70))

start.time <- Sys.time()
rslts4_u <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop, waifw = W,
  IC_type = "std", max_t = max(times_long)+100, params = paras, beep_flag = TRUE,
  adjust_beta_flag = TRUE, plot_flag = TRUE, plot_title = "realistic structure, POLYMOD (adjust beta)")
Sys.time() - start.time

start.time <- Sys.time()
rslts4_n <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop, waifw = W,
  IC_type = "std", max_t = max(times_long)+250, params = paras, beep_flag = TRUE,
  adjust_beta_flag = FALSE, plot_flag = TRUE, plot_title = "realistic structure, POLYMOD (adjust beta)")
Sys.time() - start.time

start.time <- Sys.time()
rslts4_v <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop, waifw = W,
  vax_change_times = c(0, 540, 545), vax_rates = c(0, 0.96, 0.4),
  IC_type = "std", max_t = max(times_long)+250, params = paras, beep_flag = TRUE,
  adjust_beta_flag = FALSE, plot_flag = TRUE, plot_title = "realistic structure, POLYMOD")
Sys.time() - start.time

# start.time <- Sys.time()
# rslts4_vu <- run_ode(
#   age_classes = age_classes, mort = rep(mort, length(age_classes)),
#   fert = rep(fert, length(age_classes)), start_pop = start_pop, waifw = W,
#   vax_change_times = c(0, 540, 545), vax_rates = vax_rates,
#   IC_type = "std", max_t = max(times_long)+250, params = paras, beep_flag = TRUE,
#   adjust_beta_flag = TRUE, plot_flag = TRUE, plot_title = "realistic structure, POLYMOD (adjust beta)")
# Sys.time() - start.time

# plot with and without age structure together
bind_rows(rslts4_v %>% mutate(adjust_beta = FALSE),
          rslts4_vu %>% mutate(adjust_beta = TRUE)) %>%
# rslts4_v %>% mutate(adjust_beta = FALSE) %>%
  filter(time > 500, variable %in% c("I", "BH", "BHs", "BHi")) %>%
  summarize(value = sum(value), .by = c("variable", "time", "adjust_beta")) %>%
  ggplot(aes(x = time, y = value, color = adjust_beta)) + 
  geom_line() + 
  facet_wrap(vars(variable), scales = "free", ncol = 1)

timepoints = 540:550#c(540, 542, 543.2, 545, 550)

p0 = rslts4_v %>% mutate(age = as.double(age)) %>%
  filter(time > 535, variable %in% c("S", "I")) %>%
  summarize(value = sum(value), .by = c("variable", "time")) %>%
  bind_rows(rslts4_v %>% filter(time > 535, variable == "I") %>%
              mutate(bin_width = bin_width[as.factor(age)]) %>%
              summarize(value = sum(value*(age-bin_width/2))/sum(value), .by = c("time")) %>%
              mutate(variable = "mean_age")) %>%
  ggplot() + 
  geom_rect(data = data.frame(xmin = 540, xmax = 545, ymin = -Inf, ymax = Inf), 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'gray', alpha = 0.4) + 
  geom_vline(xintercept = timepoints, linetype = "dashed", color = 'gray') +
  geom_line(aes(x = time, y = value)) + 
  facet_grid(rows = vars(variable), scales = "free") + 
  scale_x_continuous(expand = c(0,0)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank())


p1 = rslts4_v %>% mutate(age = as.double(age)) %>%
  filter(time > 535, variable %in% c("BH", "BHs", "BHi")) %>%
  summarize(value = sum(value), .by = c("variable", "time")) %>%
  ggplot() + 
  geom_hline(yintercept = 1) + 
  geom_vline(xintercept = timepoints, linetype = "dashed", color = 'gray') +
  geom_line(aes(x = time, y = value, color = variable)) + 
  geom_text(data = data.frame(y = c(Inf, -Inf), x = c(Inf, Inf), vjust = c(1, 0),
                              lab = c("age-structured transmission rate\nlower than homogeneous", "age-structured transmission rate\nhighter than homogeneous")), 
            aes(x = x, y = y, label = lab, vjust = vjust), hjust = 1) + 
  scale_x_continuous(expand = c(0,0)) + 
  theme_bw() + 
  theme(legend.position = "right", 
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank())


p2 = rslts4_v %>% filter(time %in% timepoints, variable %in% c("S", "I"), age <= 40) %>%
  mutate(bin_width = bin_width[as.factor(age)],
         value = value/bin_width) %>%
  ggplot(aes(x = age, y = value)) +
  geom_hline(data = rslts4_v %>% filter(time %in% timepoints, variable %in% c("S", "I")) %>%
               summarize(value = sum(value)/n(), .by = c("variable", "time")),
             aes(yintercept = value), color = "red") +
  geom_point(size = 1, shape = 1) +
  facet_grid(rows = vars(variable), cols = vars(time), scales = "free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

p22 = rslts4_v %>% filter(time %in% timepoints, variable %in% c("S", "I"), age <= 10) %>%
  mutate(bin_width = bin_width[as.factor(age)], 
         value = value/bin_width) %>%
  dcast(time + age ~ variable, value.var = "value") %>%
  ggplot(aes(x = S, y = I, color = as.factor(age))) + 
  geom_vline(data = rslts4_v %>% filter(time %in% timepoints, variable %in% c("S")) %>%
               summarize(value = sum(value)/n(), .by = c("variable", "time")),
             aes(xintercept = value), color = "black", linetype = "dotted") +
  geom_hline(data = rslts4_v %>% filter(time %in% timepoints, variable %in% c("I")) %>%
               summarize(value = sum(value)/n(), .by = c("variable", "time")),
             aes(yintercept = value), color = "black", linetype = "dotted") +
  geom_path(color = "gray") +
  geom_point() + 
  facet_wrap(vars(time), nrow = 1) + 
  scale_color_viridis_d(direction = -1) + 
  # scale_x_log10() +
  # scale_y_log10() +
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank())

p23 = rslts4_v %>% filter(time > 535, variable %in% c("S", "I"), age <= 30) %>%
  mutate(bin_width = bin_width[as.factor(age)], 
         value = value/bin_width) %>% 
  dcast(time + age ~ variable, value.var = "value") %>%
  ggplot(aes(x = time, y = as.factor(age))) + 
  geom_tile(aes(fill = S)) + 
  geom_point(data = rslts4_v %>% filter(time %in% 535:550, variable == "I", age <= 30) %>%
               mutate(bin_width = bin_width[as.factor(age)], 
                      I = value/bin_width), 
             aes(size = I), color = "white") + 
  scale_fill_viridis_c() + 
  scale_size_continuous(range = c(0,3)) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  theme(legend.position = "right")

cowplot::plot_grid(p0, p1, p23, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.35, 0.3, 0.35))


rslts4_v %>% filter(time %in% timepoints, variable %in% c("S", "I"), age <= 10) %>%
  mutate(bin_width = bin_width[as.factor(age)], 
         value = value/bin_width) %>%
  dcast(time + age ~ variable, value.var = "value") %>%
  ggplot(aes(x = S, y = I, color = as.factor(age))) + 
  geom_vline(data = rslts4_v %>% filter(time %in% timepoints, variable %in% c("S")) %>%
              summarize(value = sum(value)/n(), .by = c("variable", "time")),
              aes(xintercept = value), color = "black", linetype = "dotted") +
  geom_hline(data = rslts4_v %>% filter(time %in% timepoints, variable %in% c("I")) %>%
               summarize(value = sum(value)/n(), .by = c("variable", "time")),
             aes(yintercept = value), color = "black", linetype = "dotted") +
  geom_path(color = "gray") +
  geom_point() + 
  facet_wrap(vars(time), scales = "free") + 
  scale_color_viridis_d(direction = -1) + 
  # scale_x_log10() +
  # scale_y_log10() +
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank())

rslts4_v %>% filter(time %in% timepoints, variable %in% c("S"), age <= 10) %>%
  rename(age_S = age, S = value) %>%
  select(-variable) %>%
  left_join(rslts4_v %>% filter(time %in% timepoints, variable %in% c("I"), age <= 10) %>%
              rename(age_I = age, I = value) %>%
              select(-variable)) %>%
  dplyr::mutate(row_id = 1:dplyr::n()) %>%
  mutate(age_S_indx = which(age_classes == age_S),
         age_I_indx = which(age_classes == age_I), 
         contact_rate = W[age_S_indx, age_I_indx], .by = "row_id") %>%
  
  
  