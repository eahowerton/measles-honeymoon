library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidytable)
library(epimdr2)
library(cowplot)
library(dplyr)
library(beepr)
library(fields)

source("R/age-struc-functions.R")

#### MODEL SETUP ---------------------------------------------------------------
compartments = c("S", "E", "I", "R")
times = seq(0, 150, by = 1/365)
times_long = seq(0, 300, by = 1/365)

no_vax = approxfun(times_long, rep(0, length(times_long)))

paras = c(mu = 1/50, N = 1, beta0 = 1000, beta1 = 0, omega = 0,
          sigma = 365/8, gamma = 365/5, vax_pct = 0, delta = 0, p = 0)

# R0
with(as.list(paras), beta0/(gamma + mu)*(sigma/(sigma + mu)))

start_pop = 1
fert = paras["mu"]
mort = paras["mu"]

vax_rate_high = 0.8
vax_rate_drop = 0.4

#### TEST WITH ONE COMPARTMENT -------------------------------------------------
age_classes1 =  c(100)
rslts1 <- run_ode(
  age_classes = age_classes1, mort = mort, fert = fert, start_pop = start_pop, 
  IC_type = "std", max_t = 400, params = paras, compartments = compartments,
  plot_flag = TRUE, plot_title = "1 age class")

new_IC = rslts1 %>% filter(time == max(time), variable %in% c("S", "E", "I", "R")) %>% pull(value)
names(new_IC) = c("S", "E", "I", "R")

rslts1_v <- run_ode(
  age_classes = age_classes1, mort = mort, fert = fert, start_pop = start_pop, 
  IC_type = "std", max_t = 400, params = paras, compartments = compartments,
  vax_change_times = c(0), vax_rates = c(vax_rate_high),
  plot_flag = TRUE, plot_title = "1 age class")

new_IC_highvax = rslts1_v %>% filter(time == max(time), variable %in% c("S", "E", "I", "R")) %>% pull(value)
names(new_IC_highvax) = c("S", "E", "I", "R")

#### TEST WITH TWO COMPARTMENTS ------------------------------------------------
age_classes2 = c(20, 40)
rslts2 <- run_ode(
  age_classes = age_classes2, mort = rep(mort, length(age_classes2)),
  fert = rep(fert, length(age_classes2)), start_pop = start_pop, compartments = compartments,
  IC_type = "manual", IC_manual = new_IC, max_t = 50, params = paras,
  plot_flag = TRUE, plot_title = "2 age classes")

#### TEST WITH REALISTIC STRUCTURE ---------------------------------------------
# up to 5 years in months, up to 10 in years, and up to 70 in 5 years 
# (multiply by 4 to go from months to weeks)
# age_classes = c(seq(0.5, 9.5, 1/2), seq(10,70,by=4))# 6:9,
age_classes = seq(2, 70, 2)
bin_width = diff(c(0, age_classes))
names(bin_width) = age_classes

start.time <- Sys.time()
rslts3 <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop,
  IC_type = "manual", IC_manual = new_IC, max_t = 50, params = paras, beep_flag = TRUE,
  plot_flag = TRUE, plot_title = "realistic structure, constant WAIFW")
Sys.time() - start.time

#### NOW ADD VACCINATION -------------------------------------------------------
vax_change_times = c(0, 2)
vax_rates = c(vax_rate_high, vax_rate_drop)

start.time <- Sys.time()
rslts3_v <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop, 
  vax_change_times = vax_change_times, vax_rates = vax_rates,
  IC_type = "manual", IC_manual = new_IC_highvax, max_t = 50, params = paras, beep_flag = TRUE,
  plot_flag = TRUE, plot_title = "realistic structure, constant WAIFW")
Sys.time() - start.time

ggplot(data = rslts4_u %>% filter(variable == "I") %>% summarize(value = sum(value), .by = c("time"))) + 
  geom_line(aes(x = time, y = value))

#### TEST WITH POLYMOD CONTACTS ------------------------------------------------
# now add more realistic mixing with POLYMOD data
W = create_polymod_matrix(age_classes, plot_flag = TRUE, age_classes_to_label = c(2, 4, 6, 8, 10, 30, 50, 70))

start.time <- Sys.time()
rslts4_u <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop, waifw = W,
  IC_type = "manual", IC_manual = new_IC, max_t = 50, params = paras, beep_flag = TRUE,
  adjust_beta_flag = TRUE, plot_flag = TRUE, plot_title = "realistic structure, POLYMOD (adjust beta)")
Sys.time() - start.time

start.time <- Sys.time()
rslts4_n <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop, waifw = W,
  IC_type = "manual", IC_manual = new_IC, max_t = 100, params = paras, beep_flag = TRUE,
  adjust_beta_flag = FALSE, plot_flag = TRUE, plot_title = "realistic structure, POLYMOD (don't adjust beta)")
Sys.time() - start.time

start.time <- Sys.time()
rslts4_v <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop, waifw = W,
  vax_change_times = c(0, 50), vax_rates = vax_rates,
  IC_type = "manual", IC_manual = new_IC_highvax, max_t = 100, params = paras, beep_flag = TRUE,
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

# timepoints = seq(448, 458, 2)
timepoints = c(2, 4, 6)+50
max_time = 65
min_time = 45

p0 = rslts4_v %>% mutate(age = as.double(age)) %>%
  filter(time > min_time, time < max_time, variable %in% c("S", "I")) %>%
  summarize(value = sum(value), .by = c("variable", "time")) %>%
  bind_rows(rslts4_v %>% filter(time > min_time, time < max_time, variable == "I") %>%
              mutate(bin_width = bin_width[as.factor(age)]) %>%
              summarize(value = sum(value*(age-bin_width/2))/sum(value), .by = c("time")) %>%
              mutate(variable = "mean_age")) %>%
  ggplot() + 
  # geom_rect(data = data.frame(xmin = 540, xmax = 545, ymin = -Inf, ymax = Inf), 
  #           aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'gray', alpha = 0.4) + 
  # geom_vline(xintercept = timepoints, linetype = "dashed", color = 'gray') +
  geom_line(aes(x = time, y = value)) + 
  facet_grid(rows = vars(variable), scales = "free") + 
  scale_x_continuous(expand = c(0,0)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank())


p1 = rslts4_v %>% mutate(age = as.double(age)) %>%
  filter(time > min_time, time < max_time, variable %in% c("BH", "BHs", "BHi")) %>%
  summarize(value = sum(value), .by = c("variable", "time")) %>%
  ggplot() + 
  geom_hline(yintercept = 1) + 
  # geom_vline(xintercept = timepoints, linetype = "dashed", color = 'gray') +
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
# 
# p22 = rslts4_v %>% filter(time %in% timepoints, variable %in% c("S", "I"), age <= 10) %>%
#   mutate(bin_width = bin_width[as.factor(age)], 
#          value = value/bin_width) %>%
#   dcast(time + age ~ variable, value.var = "value") %>%
#   ggplot(aes(x = S, y = I, color = as.factor(age))) + 
#   geom_vline(data = rslts4_v %>% filter(time %in% timepoints, variable %in% c("S")) %>%
#                summarize(value = sum(value)/n(), .by = c("variable", "time")),
#              aes(xintercept = value), color = "black", linetype = "dotted") +
#   geom_hline(data = rslts4_v %>% filter(time %in% timepoints, variable %in% c("I")) %>%
#                summarize(value = sum(value)/n(), .by = c("variable", "time")),
#              aes(yintercept = value), color = "black", linetype = "dotted") +
#   geom_path(color = "gray") +
#   geom_point() + 
#   facet_wrap(vars(time), nrow = 1) + 
#   scale_color_viridis_d(direction = -1) + 
#   # scale_x_log10() +
#   # scale_y_log10() +
#   theme_bw() + 
#   theme(legend.position = "bottom", 
#         panel.grid.minor = element_blank(), 
#         panel.grid.major.x = element_blank())

p23 = rslts4_v %>% filter(time > min_time, time < max_time, variable %in% c("S", "I"), age <= 30) %>%
  mutate(bin_width = bin_width[as.factor(age)], 
         value = value/bin_width) %>% 
  dcast(time + age ~ variable, value.var = "value") %>%
  ggplot(aes(x = time, y = as.factor(age))) + 
  geom_tile(aes(fill = S)) + 
  geom_point(data = rslts4_v %>% filter(time %in% seq(min_time, max_time, 4), variable == "I", age <= 30) %>%
               mutate(bin_width = bin_width[as.factor(age)],
                      I = value/bin_width),
             aes(size = I), color = "white") +
  scale_fill_viridis_c() + 
  scale_size_continuous(range = c(0,3)) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  theme(legend.position = "right")

cowplot::plot_grid(p0, p1, p23, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.35, 0.3, 0.35))


rslts4_v %>% filter(time %in% timepoints, variable %in% c("S", "I")) %>%
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
  facet_wrap(vars(time), ncol = 2, scales = "free") +
  scale_color_viridis_d(direction = -1) + 
  # scale_x_log10() +
  # scale_y_log10() +
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank())

p24 = rslts4_v %>% filter(time %in% 540:550, variable %in% c("S"), age <= 10) %>%
  rename(age_S = age, S = value) %>%
  select(-variable) %>%
  left_join(rslts4_v %>% filter(time %in% timepoints, variable %in% c("I"), age <= 10) %>%
              rename(age_I = age, I = value) %>%
              select(-variable)) %>%
  dplyr::mutate(row_id = 1:dplyr::n()) %>%
  mutate(age_S_indx = which(age_classes == age_S),
         age_I_indx = which(age_classes == age_I), 
         contact_rate = W[age_S_indx, age_I_indx], .by = "row_id") %>%
  summarize(weighted_S = sum(S*contact_rate), .by = c("age_I", "I", "time")) %>%
  mutate(bin_width = bin_width[as.factor(age_I)], 
         I = I/bin_width) %>%
  ggplot(aes(x = age_I, y = I, color = weighted_S)) + 
  geom_point() + 
  facet_wrap(vars(time), nrow = 1) + 
  scale_color_viridis_c() + 
  theme_bw()
  
  
#### CALCULATE EQUILIBRIUM UNITY BETA AT DIFFERENT VACCINATION COVERAGES -------
test_v = seq(0, 0.75, 0.15)
equil_vax = vector("list", length(test_v))
for(i in 1:length(test_v)){
  print(paste0(i, "/", length(test_v)))
  equil_vax[[i]] = run_ode(
    age_classes = age_classes, mort = rep(mort, length(age_classes)),
    fert = rep(fert, length(age_classes)), start_pop = start_pop, waifw = W,
    compartments = compartments,
    vax_change_times = c(0), vax_rates = c(test_v[i]),
    IC_type = "manual", IC_manual = new_IC, max_t = ifelse(i < 4, 200, 400), params = paras, 
    adjust_beta_flag = FALSE, plot_flag = FALSE) %>%
    mutate(v = test_v[i])
}

saveRDS(equil_vax, "data/output-data/equil_vax.rds")

test_v = seq(0, 0.75, 0.15)
equil_vax_flat = vector("list", length(test_v))
for(i in 1:length(test_v)){
  print(paste0(i, "/", length(test_v)))
  equil_vax_flat[[i]] = run_ode(
    age_classes = age_classes, mort = rep(mort, length(age_classes)),
    fert = rep(fert, length(age_classes)), start_pop = start_pop,
    compartments = compartments,
    vax_change_times = c(0), vax_rates = c(test_v[i]),
    IC_type = "manual", IC_manual = new_IC, max_t = ifelse(i < 4, 200, 350), params = paras, 
    adjust_beta_flag = FALSE, plot_flag = FALSE) %>%
    mutate(v = test_v[i])
}

equil_vax_full = bind_rows(equil_vax)
equil_vax_flat_full = bind_rows(equil_vax_flat)

# beta hats
bind_rows(
  equil_vax_full %>%
    mutate(max_time = max(time), .by = c("v")) %>%
    filter(time == max_time, variable %in% c("BH", "BHi", "BHs")) %>%
    mutate(w = "POLYMOD"), 
  equil_vax_flat_full %>%
    mutate(max_time = max(time), .by = c("v")) %>%
    filter(time == max_time, variable %in% c("BH", "BHi", "BHs")) %>%
    mutate(w = "flat")
) %>%
  ggplot(aes(x = v, y = value, color = variable)) + 
  geom_hline(yintercept = 1, color = "black") + 
  geom_text(data = data.frame(y = c(Inf, -Inf), x = c(Inf, Inf), vjust = c(1, 0),
                              lab = c("age-structured transmission rate\nlower than homogeneous\n(slowing down with age structure)", "")), 
                                      #"(speeding up with age structure)\nage-structured transmission rate\nhighter than homogeneous")), 
            aes(x = x, y = y, label = lab, vjust = vjust), hjust = 1, color = "black") + 
  geom_point() + 
  geom_line(aes(linetype = w)) + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  theme_bw()

# equilibrium infections
bind_rows(
  equil_vax_full %>%
  mutate(max_time = max(time), .by = c("v")) %>%
  filter(time == max_time, variable %in% c("I")) %>% 
    mutate(w = "POLYMOD"), 
  equil_vax_flat_full %>%
    mutate(max_time = max(time), .by = c("v")) %>%
    filter(time == max_time, variable %in% c("I")) %>% 
    mutate(w = "flat")) %>%
  summarize(value = sum(value), .by = c("v", "w")) %>%
  ggplot(aes(x = v, y = value)) + 
  geom_point() + 
  geom_line(aes(linetype = w)) + 
  labs(x = "vaccination coverage", y = "equilibrium I") +
  scale_linetype_manual(values = c("dashed", "solid")) + 
  theme_bw()

# equilibrium mean age
bind_rows(
  equil_vax_full %>% 
  mutate(max_time = max(time), .by = c("v")) %>%
  filter(time == max_time, variable %in% c("I")) %>%
    mutate(w = "POLYMOD"), 
  equil_vax_flat_full %>% 
    mutate(max_time = max(time), .by = c("v")) %>%
    filter(time == max_time, variable %in% c("I")) %>%
    mutate(w = "flat")) %>%
  mutate(bin_width = bin_width[as.factor(age)]) %>%
  summarize(mean_age = sum(value*(age-bin_width/2))/sum(value), .by = c("time", "v", "w")) %>%
  ggplot(aes(x = v, y = mean_age)) +
  geom_point() + 
  geom_line(aes(linetype = w)) + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  theme_bw()

# age distributions of S and I
equil_vax_full %>%
  mutate(max_time = max(time), .by = c("v")) %>%
  filter(time == max_time, variable %in% c("S", "I")) %>%
  ggplot() + 
  geom_line(aes(x = age, y = value)) + 
  geom_hline(data = equil_vax_full %>%
               mutate(max_time = max(time), .by = c("v")) %>%
               filter(time == max_time, variable %in% c("S", "I")) %>%
              summarize(mean = mean(value), .by = c("v", "variable")), 
            aes(yintercept = mean), color = "red") +
  facet_grid(rows = vars(variable), cols = vars(v), scales = "free") + 
  theme_bw()


#### NOW START AT VAX EQUILIBRIUM AND RELEASE ----------------------------------
new_IC_polymod_vax = equil_vax_full %>% 
  filter(time == max(time), variable %in% c("S", "E", "I", "R"), v == 0.6) %>% 
  arrange(age) %>%
  mutate(name = paste0(variable, "_", age))
new_IC_polymod_vax_name = new_IC_polymod_vax %>% pull(name)
new_IC_polymod_vax = new_IC_polymod_vax %>% pull(value)
names(new_IC_polymod_vax) = new_IC_polymod_vax_name

new_IC_polymod_vax_flat = equil_vax_flat_full %>% 
  filter(time == max(time), variable %in% c("S", "E", "I", "R"), v == 0.6) %>% 
  arrange(age) %>%
  mutate(name = paste0(variable, "_", age))
new_IC_polymod_vax_flat_name = new_IC_polymod_vax_flat %>% pull(name)
new_IC_polymod_vax_flat = new_IC_polymod_vax_flat %>% pull(value)
names(new_IC_polymod_vax_flat) = new_IC_polymod_vax_flat_name


equil_vax_p6 = run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop, waifw = W,
  compartments = compartments,
  vax_change_times = c(0), vax_rates = c(0.6),
  IC_type = "manual", IC_manual = new_IC_polymod_vax, max_t = 100, params = paras, 
  adjust_beta_flag = FALSE, plot_flag = FALSE)

new_IC_polymod_vax = equil_vax_p6 %>% 
  filter(time == max(time), variable %in% c("S", "E", "I", "R")) %>% 
  arrange(age) %>%
  mutate(name = paste0(variable, "_", age))
new_IC_polymod_vax_name = new_IC_polymod_vax %>% pull(name)
new_IC_polymod_vax = new_IC_polymod_vax %>% pull(value)
names(new_IC_polymod_vax) = new_IC_polymod_vax_name


equil_vax_release = run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop,
  compartments = compartments,
  vax_change_times = c(0, 2), vax_rates = c(0.6, 0.1), waifw = W,
  IC_type = "manual", IC_manual = new_IC_polymod_vax, max_t = 25, params = paras, 
  adjust_beta_flag = FALSE, plot_flag = FALSE) %>%
  mutate(v = test_v[i])

equil_vax_flat_release = run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop,
  compartments = compartments,
  vax_change_times = c(0, 2), vax_rates = c(0.6, 0.1),
  IC_type = "manual", IC_manual = new_IC_polymod_vax_flat, max_t = 25, params = paras, 
  adjust_beta_flag = FALSE, plot_flag = FALSE) %>%
  mutate(v = test_v[i])


# beta hats
p1 = bind_rows(
  equil_vax_release %>%
    filter(variable %in% c("BH", "BHi", "BHs")) %>%
    mutate(w = "POLYMOD"),
  equil_vax_flat_release %>%
    filter(variable %in% c("BH", "BHi", "BHs")) %>%
    mutate(w = "flat")
) %>%
  ggplot(aes(x = time, y = value, color = variable)) + 
  geom_vline(xintercept = 2, color = "black") +
  geom_hline(yintercept = 1, color = "black") + 
  geom_text(data = data.frame(y = c(Inf, -Inf), x = c(Inf, Inf), vjust = c(1, 0),
                              lab = c("age-structured transmission rate\nlower than homogeneous\n(slowing down with age structure)", "")), 
            #"(speeding up with age structure)\nage-structured transmission rate\nhighter than homogeneous")), 
            aes(x = x, y = y, label = lab, vjust = vjust), hjust = 1, color = "black") + 
  # geom_point() + 
  geom_line(aes(linetype = w)) + 
  labs(y = "unity beta") +
  scale_linetype_manual(values = c("dashed", "solid")) + 
  theme_bw()
# infections
p2 = bind_rows(
  equil_vax_release %>%
    filter(variable %in% c("I")) %>%
    mutate(w = "POLYMOD"),
  equil_vax_flat_release %>%
    filter(variable %in% c("I")) %>%
    mutate(w = "flat")
) %>%
  summarize(value = sum(value), .by = c("time", "w")) %>%
  ggplot(aes(x = time, y = value)) + 
  geom_vline(xintercept = 2, color = "black") +
  geom_line(aes(linetype = w)) + 
  labs(y = "infections") + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  theme_bw()
plot_grid(p1, p2, ncol = 1, align = "v", axis ="lr")

#### TRY WITH JESS'S MATRICES --------------------------------------------------
xx = findStableStruct(age_classes, fert = rep(fert, length(age_classes)), 
                      mort = rep(mort, length(age_classes)))

# measles like?
jcm_results = vector("list", length(waifw))
paras_jcm = c(mu = 1/50, N = 1, beta0 = 365, beta1 = 0, omega = 0,
                      sigma = 365/1, gamma = 365/14, vax_pct = 0, delta = 0, p = 0)
paras_jcm["N"] = 500000

# R0
with(as.list(paras_jcm), beta0/(gamma + mu)*(sigma/(sigma + mu)))
# flatW = scale.Waifw(R0=14, DFE=xx$stable.age*1e6, waifw=matrix(0.0002,length(age.classes),length(age.classes)))

# re-scale Jess's WAIFW matrices to work with my parameters
waifw2 = waifw
waifw2[[5]] = unname(waifw2[[5]])
waifw2 = lapply(waifw2, function(i){i/mean(i)})

lapply(waifw2, melt) %>%
  bind_rows(.id = "waifw_id") %>%
  mutate(value_scaled = value/max(value), .by = c("waifw_id")) %>%
  ggplot(aes(x = Var1, y = Var2, fill = value_scaled)) + 
  geom_tile() + 
  facet_wrap(vars(waifw_id)) +
  scale_fill_viridis_c() + 
  scale_x_continuous(expand = c(0,0),
                     breaks = which(age_classes %in% c(2, 4, 6, 8, 10, 30, 50, 70)),
                     labels = c(2, 4, 6, 8, 10, 30, 50, 70),
                     name = "age (years)") +
  scale_y_continuous(expand = c(0,0),
                     breaks = which(age_classes %in% c(2, 4, 6, 8, 10, 30, 50, 70)),
                     labels = c(2, 4, 6, 8, 10, 30, 50, 70),
                     name = "age (years)")

# do the matching on age-specific incidence 
scalars = match_on_age_or_inc(age_classes = age_classes, mort =  rep(mort, length(age_classes)), 
                    fert = rep(fert, length(age_classes)), 
                    start_pop = paras_jcm["N"], compartments = compartments, 
                    params = paras_jcm, waifws = waifw2[-1], plot_flag = TRUE)

chosen_scalars = scalars %>% filter(variable == "tot_I_diff", best_value == TRUE)



for(i in 1:length(waifw)){
  print(paste0(i, "/", length(waifw)))
  paras_tmp = paras_jcm
  if(i > 1){
    paras_tmp["beta0"] = paras_tmp["beta0"] * chosen_scalars %>% filter(waifw_id == i) %>% pull(scalar)
  }
  jcm_results[[i]] = run_ode(
    age_classes = age_classes, mort = rep(mort, length(age_classes)),
    fert = rep(fert, length(age_classes)), start_pop = paras_jcm["N"],
    compartments = compartments, w = waifw2[[i]],
    IC_type = "std", max_t = 100, params = paras_tmp, #, IC_manual = new_IC
    adjust_beta_flag = FALSE, plot_flag = FALSE) #%>%
    # mutate(v = test_v[i])
}

jcm_results_long = bind_rows(jcm_results, .id = "waifw_id")

jcm_results_long %>%
  filter(variable %in% c("I")) %>%
  summarize(value = sum(value), .by = c("time", "variable", "waifw_id")) %>%
  ggplot(aes(x = time, y = value, color = as.factor(waifw_id))) + 
  geom_line() + 
  facet_wrap(vars(variable), scales = "free")

jcm_results_long %>%
  filter(variable %in% c("BH", "BHs", "BHi"), time >10) %>%
  summarize(value = sum(value), .by = c("time", "variable", "waifw_id")) %>%
  ggplot(aes(x = time, y = value, color = variable)) + 
  geom_hline(yintercept = 1) + 
  geom_line() + 
  facet_wrap(vars(waifw_id)) + 
  theme_bw()

jcm_results_long %>%
  filter(variable %in% c("I"), time == 1) %>%
  ggplot(aes(x = age, y = value))+ 
  geom_line() + 
  facet_wrap(vars(waifw_id), scales = "free")

