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
library(patchwork)

source("R/age-struc-functions.R")

#### MODEL SETUP ---------------------------------------------------------------
compartments = c("S", "E", "I", "R")
times = seq(0, 150, by = 1/365)
times_long = seq(0, 300, by = 1/365)

no_vax = approxfun(times_long, rep(0, length(times_long)))

paras = c(mu = 1/50, N = 1, beta0 = 1000, beta1 = 0, omega = 0, delta = 0,
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
  age_classes = age_classes1, mort = rep(paras["mu"], length(age_classes1)), 
  fert = rep(paras["mu"], length(age_classes1)), start_pop = start_pop, 
  IC_type = "std", max_t = 400, params = paras, compartments = compartments,
  plot_flag = TRUE, plot_title = "1 age class")

new_IC = rslts1 %>% filter(time == max(time), variable %in% c("S", "E", "I", "R")) %>% pull(value)
names(new_IC) = c("S", "E", "I", "R")

rslts1_v <- run_ode(
  age_classes = age_classes1, mort = rep(paras["mu"], length(age_classes1)), 
  fert = rep(paras["mu"], length(age_classes1)), start_pop = start_pop, 
  IC_type = "std", max_t = 400, params = paras, compartments = compartments,
  vax_change_times = c(0), vax_rates = c(vax_rate_high),
  plot_flag = TRUE, plot_title = "1 age class")

new_IC_highvax = rslts1_v %>% filter(time == max(time), variable %in% c("S", "E", "I", "R")) %>% pull(value)
names(new_IC_highvax) = c("S", "E", "I", "R")

#### TEST WITH TWO COMPARTMENTS ------------------------------------------------
age_classes2 = c(20, 40)
rslts2 <- run_ode(
  age_classes = age_classes2, mort = rep(paras["mu"], length(age_classes2)),
  fert = rep(paras["mu"], length(age_classes2)), start_pop = start_pop, compartments = compartments,
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
  age_classes = age_classes, mort = rep(paras["mu"], length(age_classes)),
  compartments = compartments,
  fert = rep(paras["mu"], length(age_classes)), start_pop = start_pop,
  IC_type = "manual", IC_manual = new_IC, max_t = 50, params = paras, beep_flag = TRUE,
  plot_flag = TRUE, plot_title = "realistic structure, constant WAIFW")
Sys.time() - start.time

#### NOW ADD VACCINATION -------------------------------------------------------
vax_change_times = c(0, 2)
vax_rates = c(vax_rate_high, vax_rate_drop)

start.time <- Sys.time()
rslts3_v <- run_ode(
  age_classes = age_classes, mort = rep(paras["mu"], length(age_classes)),
  fert = rep(paras["mu"], length(age_classes)), start_pop = start_pop, 
  vax_change_times = vax_change_times, vax_rates = vax_rates,
  IC_type = "manual", IC_manual = new_IC_highvax, max_t = 50, params = paras, beep_flag = TRUE,
  plot_flag = TRUE, plot_title = "realistic structure, constant WAIFW")
Sys.time() - start.time

ggplot(data = rslts4_u %>% filter(variable == "I") %>% summarize(value = sum(value), .by = c("time")) %>% mutate(convert_beta = TRUE) %>%
         # bind_rows(rslts4_n %>% filter(variable == "I") %>% summarize(value = sum(value), .by = c("time")) %>% mutate(convert_beta = FALSE)) %>%
         bind_rows(rslts1 %>% filter(variable == "I") %>% mutate(convert_beta = "baseline")) %>%
         filter(time < 50)
         ) + 
  geom_line(aes(x = time, y = value, color = convert_beta))

#### TEST WITH POLYMOD CONTACTS ------------------------------------------------
# now add more realistic mixing with POLYMOD data
W = create_polymod_matrix(age_classes, plot_flag = TRUE, age_classes_to_label = c(2, 4, 6, 8, 10, 30, 50, 70))

start.time <- Sys.time()
rslts4_u <- run_ode(
  age_classes = age_classes, mort = rep(paras["mu"], length(age_classes)),
  compartments = compartments, 
  fert = rep(paras["mu"], length(age_classes)), start_pop = start_pop, waifw = W,
  IC_type = "manual", IC_manual = new_IC, max_t = 200, params = paras, beep_flag = TRUE,
  adjust_beta_flag = TRUE, plot_flag = TRUE, plot_title = "realistic structure, POLYMOD (adjust beta)")
Sys.time() - start.time

start.time <- Sys.time()
rslts4_n <- run_ode(
  age_classes = age_classes_jcm, mort = rep(paras["mu"], length(age_classes_jcm)),
  compartments = compartments, 
  fert = rep(paras["mu"], length(age_classes_jcm)), start_pop = start_pop, waifw = W,
  IC_type = "manual", IC_manual = new_IC, max_t = 100, params = paras, beep_flag = TRUE,
  adjust_beta_flag = FALSE, plot_flag = TRUE, plot_title = "realistic structure, POLYMOD (don't adjust beta)")
Sys.time() - start.time

start.time <- Sys.time()
rslts4_v <- run_ode(
  age_classes = age_classes, mort = rep(paras["mu"], length(age_classes)),
  fert = rep(paras["mu"], length(age_classes)), start_pop = start_pop, waifw = W,
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
age_classes_jcm = c(seq(4, 180, 4),seq(240,1200,by=60))/12
bin_width_jcm = diff(c(0, age_classes_jcm))

mort <- c(rep(0,length(age_classes_jcm)-1),1)
fert <-  c(rep(0,length(age_classes_jcm)-1),1)
xx = findStableStruct(age_classes_jcm, fert = fert,
                      mort = mort)

# xx = findStableStruct(age_classes_jcm, fert = rep(fert, length(age_classes_jcm)), 
#                       mort = rep(mort, length(age_classes_jcm)))

## series from flat to peaky to diagonal to polymod
waifw2 = get_waifws(age_classes_jcm)

# measles like?
jcm_results = vector("list", length(waifw))
paras_jcm = c(mu = 1/50, N = 1, beta0 = 365, beta1 = 0, omega = 0,
                      sigma = 365/1, gamma = 365/14, vax_pct = 0, delta = 0, p = 0)
paras_jcm["N"] = 500000
paras_jcm["delta"] = 1e-4

# R0
with(as.list(paras_jcm), beta0/(gamma + mu)*(sigma/(sigma + mu)))

waifw_labs = c("flat", "peak age 5", "peak age 10", "diagonal", "POLYMOD")
names(waifw_labs) = 1:5

lapply(waifw2, melt) %>%
  bind_rows(.id = "waifw_id") %>%
  mutate(value_scaled = value/max(value), .by = c("waifw_id")) %>%
  ggplot(aes(x = Var1, y = Var2, fill = value_scaled)) + 
  geom_tile() + 
  facet_wrap(vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) +
  scale_fill_viridis_c() + 
  scale_x_continuous(expand = c(0,0),
                     # breaks = which(age_classes_jcm %in% c(2, 4, 6, 8, 10, 30, 50, 70)),
                     # labels = c(2, 4, 6, 8, 10, 30, 50, 70),
                     name = "age (years)") +
  scale_y_continuous(expand = c(0,0),
                     # breaks = which(age_classes_jcm %in% c(2, 4, 6, 8, 10, 30, 50, 70)),
                     # labels = c(2, 4, 6, 8, 10, 30, 50, 70),
                     name = "age (years)")

# do the matching on age-specific incidence
scalars_v = vector("list", length(test_v))
for(i in 1:length(test_v)){
  print(paste0("V = ", test_v[i]))
  scalars_v[[i]] = match_on_age_or_inc(age_classes = age_classes_jcm, mort = mort, 
                                       fert = fert, start_pop = paras_jcm["N"], 
                                       compartments = compartments, 
                                       vax_pct = test_v[i], max_t = 200, 
                                       params = paras_jcm, waifws = waifw2[-1])
}

scalars_v_full =  bind_rows(scalars_v, .id = "vax_id") %>%
  mutate(vax_id = as.integer(vax_id)) %>%
  left_join(data.frame(vax_id = 1:length(test_v), 
                       v = test_v))

saveRDS(scalars_v_full, "data/output-data/scalars_by_v.rds")

ggplot(data = scalars_v_full, aes(x = v, y = best_scalar))  +
  geom_line() +
  geom_point() + 
  facet_wrap(vars(waifw_id))

# ggplot(data = scalars_v_full %>% filter(variable == "tot_I_diff"), aes(x = scalar, y = value)) + 
#   geom_line() + 
#   geom_point(aes(color = as.factor(best_value)), fill = "white") + 
#   facet_grid(cols = vars(waifw_id), rows = vars(variable, v), scales = "free") +
#   scale_color_manual(values = c("black", "red")) + 
#   theme_bw()

chosen_scalars = scalars_v_full 

# to equilibrium
jcm_all_v = vector("list", length(test_v))

for(j in 1:length(test_v)){
  print(paste0(j, "/", length(test_v)))
  tmp = vector("list", length(waifw))
  for(i in 1:length(waifw)){
    paras_tmp = paras_jcm
    if(i > 1){
      s = chosen_scalars %>% filter(waifw_id == i, v == test_v[j]) %>% pull(best_scalar)
      if(is.na(s)){print(paste0("NA SCALAR, SKIPPING (i = ", i)); next}
      paras_tmp["beta0"] = paras_tmp["beta0"] * s
    }
    IC = setup_IC(start_pop = paras_tmp["N"], age_classes_jcm, compartments, fert = fert, 
                  mort = mort, IC_type = "stable-age") 
    Fmat <- buildFMatrix(age.classes = age_classes_jcm, fert = fert, ncompartments = length(compartments), time.step = 1)
    tmp_steady = runsteady(y = IC, times = c(0, 500), func = sir_age_structured, parms = paras_tmp, # runsteady arguments
                       compartments = compartments, age_classes = age_classes_jcm, mort = mort, fert = fert,
                       vax_change_times = c(0), vax_rates = c(test_v[j]), waifw = waifw2[[i]],
                       Fmat = Fmat, adjust_beta_flag = FALSE, print_warnings_flag = FALSE
      )
    tmp[[i]] = data.frame(variable = names(tmp_steady$y), 
                          value = tmp_steady$y, 
                          time = attributes(tmp_steady)$time,
                          v = test_v[j], waifw_id = i, s = ifelse(i == 1, 1, s))
  }
  jcm_all_v[[j]] = bind_rows(tmp)
}
beep()

# run for a few time steps to get beta hat values
jcm_bh = vector("list", length(test_v))

for(j in 1:length(test_v)){
  print(paste0(j, "/", length(test_v)))
  tmp = vector("list", length(waifw))
  for(i in 1:length(waifw)){
    paras_tmp = paras_jcm
    if(i > 1){
      s = chosen_scalars %>% filter(waifw_id == i, v == test_v[j]) %>% pull(best_scalar)
      if(is.na(s)){print(paste0("NA SCALAR, SKIPPING (i = ", i)); next}
      paras_tmp["beta0"] = paras_tmp["beta0"] * s
    }
    IC_manual = jcm_all_v_long %>% filter(v == test_v[j], waifw_id == i, !(variable %in% c("BH", "BHs", "BHi", "BHb")))
    names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
    IC_manual = IC_manual$value
    names(IC_manual) = names_IC
    tmp[[i]] = run_ode(
      age_classes = age_classes_jcm, mort = mort,
      fert = fert, start_pop = paras_jcm["N"],
      compartments = compartments, dt = 1/365,
      vax_change_times = c(0), vax_rates = c(test_v[j]), waifw = waifw2[[i]],
      IC_type = "manual", IC_manual = IC_manual, max_t = 2, params = paras_tmp, #, IC_manual = new_IC
      adjust_beta_flag = FALSE, plot_flag = FALSE) %>%
      mutate(waifw_id = i, v = test_v[j], s = ifelse(i == 1, 1, s))
  }
  jcm_bh[[j]] = bind_rows(tmp)
}
beep()

bind_rows(jcm_bh) %>%
  filter(substr(variable,1,1) == "B", time == max(time)) %>%
  ggplot(aes(x = v, y = log(value), color = variable)) + 
  geom_hline(yintercept = 0) +
  geom_line() + 
  facet_wrap(vars(waifw_id))

bind_rows(jcm_bh) %>%
  filter(substr(variable,1,1) == "B", time == max(time)) %>%
  ggplot(aes(x = v, y = log(value), color = as.factor(waifw_id))) + 
  geom_hline(yintercept = 0) +
  geom_line() + 
  scale_color_brewer(palette = "Set1") +
  facet_wrap(vars(variable))

bind_rows(jcm_bh) %>%
  filter(substr(variable,1,1) == "S", time == max(time)) %>%
  summarize(value = sum(value), .by = c('time', "variable", "v", "waifw_id")) %>%
  ggplot(aes(x = v, y = value, color = variable)) + 
  geom_hline(yintercept = 0) +
  geom_line() + 
  facet_wrap(vars(waifw_id))

ggplot(data = chosen_scalars, aes(x = v, y = diff, color = as.factor(waifw_id))) + 
  geom_line()


bind_rows(jcm_bh) %>%
  filter(substr(variable,1,1) == "I", time == max(time)) %>%
  summarize(value = sum(value), .by = c('time', "variable", "v", "waifw_id")) %>%
  ggplot(aes(x = v, y = value, color = as.factor(waifw_id))) + 
  geom_line() 

bind_rows(jcm_bh) %>%
  filter(substr(variable,1,1) == "S") %>%
  summarize(value = sum(value), .by = c('time', "variable", "v", "waifw_id")) %>%
  ggplot(aes(x = time, y = value, color = as.factor(waifw_id))) + 
  geom_line() +
  facet_wrap(vars(v))
  
  mutate(age = ifelse(substr(variable, 1, 1) == "B", NA, as.double(substr(variable, 3, nchar(variable)))), 
         variable = ifelse(substr(variable, 1, 1) == "B", substr(variable, 1, 3), substr(variable, 1, 1)))



bind_rows(jcm_all_v) %>%
  mutate(age = ifelse(substr(variable, 1, 1) == "B", NA, as.double(substr(variable, 3, nchar(variable)))), 
         variable = ifelse(substr(variable, 1, 1) == "B", substr(variable, 1, 3), substr(variable, 1, 1)), 
         fert = "const") %>%
  bind_rows(bind_rows(jcm_all_v_lastd) %>%
              mutate(age = ifelse(substr(variable, 1, 1) == "B", NA, as.double(substr(variable, 3, nchar(variable)))), 
                     variable = ifelse(substr(variable, 1, 1) == "B", substr(variable, 1, 3), substr(variable, 1, 1)), 
                     fert = "struc")) %>%
  filter(variable %in% c("S", "E", "I", "R")) %>%
  mutate(age = round(age, 3)) %>%
  left_join(data.frame(age = round(age_classes_jcm,3), 
                       bin_width = bin_width_jcm)) %>%
  ggplot(aes(x = age, y = value/bin_width, color = as.factor(fert))) + 
  geom_line() + 
  facet_grid(cols = vars(v), rows = vars(variable), scales = "free")


jcm_all_v_long = bind_rows(jcm_all_v) %>%
  mutate(age = ifelse(substr(variable, 1, 1) == "B", NA, as.double(substr(variable, 3, nchar(variable)))), 
         variable = ifelse(substr(variable, 1, 1) == "B", substr(variable, 1, 3), substr(variable, 1, 1)))

jcm_all_v_long %>%
  filter(variable %in% c("I")) %>%
  summarize(value = sum(value), .by = c("variable", "waifw_id", "v")) %>%
  ggplot(aes(x = v, y = value, color = as.factor(waifw_id))) + 
  geom_point() + 
  geom_line() + 
  facet_grid(cols = vars(variable), scales = "free") + 
  scale_color_brewer(palette = "Set1", labels = waifw_labs)  + 
  theme_bw()

p1 = jcm_all_v_long %>% filter(variable == "I") %>%
  mutate(age = round(age,3)) %>%
  left_join(data.frame(age = round(age_classes_jcm,3), 
                       bin_width = bin_width_jcm)) %>%
  mutate(midpoint = age - bin_width/2) %>%
  summarize(mean_age = sum(value*midpoint*bin_width)/sum(value*bin_width), 
            var = sum((midpoint - sum(value*midpoint*bin_width)/sum(value*bin_width))^2*value*bin_width)/sum(value*bin_width),
            .by = c("v", "waifw_id")) %>%
  mutate(cv = sqrt(var)/mean_age) %>%
  melt(c("v", "waifw_id")) %>%
  filter(variable != "var") %>%
  ggplot(aes(x = v, y = value, color = as.factor(waifw_id))) + 
  geom_line(size = 0.8) + 
  facet_grid(rows = vars(variable), scales = "free") +
  scale_color_brewer(palette = "Set1", labels = waifw_labs)  + 
  scale_y_log10() +
  theme_bw() 
p1

p2 = jcm_all_v_long %>%
  filter(variable %in% c("BH", "BHs", "BHi", "BHb")) %>%
  summarize(value = sum(value)/time, .by = c("variable", "waifw_id", "v")) %>%
  ggplot(aes(x = v, y = log(value), color = variable, size = variable)) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = c(0.3, 0.75), linetype = "dotted") + 
  geom_line() + 
  geom_text(data = data.frame(y = c(Inf, -Inf), x = c(0.5,0.5), vjust = c(1, 0),
                              waifw_id = 1,
                              lab = c("\nage-structured beta\nlower than homogeneous\n(slowing down)", "(speeding up)\nage-structured beta\nhigher than homogeneous\n")), 
            #"(speeding up with age structure)\nage-structured transmission rate\nhighter than homogeneous")), 
            aes(x = x, y = y, label = lab, vjust = vjust),color = "black", size = 3) + 
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_color_manual(values = c("black", "gray", "red", "blue")) + 
  scale_size_manual(values = c(1, 0.8, 0.8, 0.8)) +
  scale_y_continuous(limits = c(-3.5,3.5), name = "log(beta hat)") +
  theme_bw() + 
  theme(panel.grid = element_blank())

# now what about waifw 2 vs. waifw 3
p3 = jcm_all_v_long %>%
  filter(v %in% c(0.3, 0.75, 0.92), variable %in% c("S", "I")) %>%
  mutate(age = round(age,3)) %>%
  left_join(data.frame(age = round(age_classes_jcm,3), 
                       bin_width = bin_width_jcm)) %>%
  ggplot(aes(x = age, y = value/bin_width, color = as.factor(waifw_id), alpha = as.factor(v))) +
  geom_line(data = jcm_all_v_long %>%
               filter(v %in% c(0.3, 0.75, 0.92), variable %in% c("S", "I")) %>%
               mutate(age = round(age,3)) %>%
               left_join(data.frame(age = round(age_classes_jcm,3), 
                                    bin_width = bin_width_jcm)) %>%
              mutate(tot_val = sum(value), 
                     tot_binwidth = sum(bin_width), .by = c("variable", "v")) %>%
              mutate(const_value = tot_val * bin_width/tot_binwidth),
             aes(x = age, y = const_value/bin_width, alpha = as.factor(v)), color = "black") +
  geom_line(size = 0.8) + 
  facet_grid(cols = vars(waifw_id), rows = vars(variable), scales = "free", 
             labeller = labeller(waifw_id = waifw_labs)) + 
  scale_alpha_manual(values = c(0.4, 0.7, 1)) +
  scale_color_brewer(palette = "Set1", labels = waifw_labs)  + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank())
p3

cowplot::plot_grid(p1, p2, p3, align = "v", axis = "lr", ncol = 1, rel_heights = c(0.2, 0.15, 0.2))

jcm_all_v_long %>%
  # filter(age < 100) %>%
  summarize(value = sum(value), .by = c("variable", "waifw_id", "v")) %>%
  filter(variable %in% c("S", "I")) %>%
  dcast(waifw_id + v ~ variable) %>%
  ggplot(aes(x = S, y = I)) + 
  geom_path(aes(group = as.factor(waifw_id)), color = "gray") +
  geom_point(aes(color = as.factor(v))) + 
  facet_wrap(vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) +
  scale_color_viridis_d() + 
  theme_bw()

jcm_all_v_long %>%
  filter(age < 100) %>%
  summarize(value = sum(value), .by = c("variable", "waifw_id", "v")) %>%
  filter(variable == "S") %>%
  ggplot(aes(x = v, y = value/paras_jcm["N"])) + 
  geom_line(aes(color = as.factor(waifw_id))) +
  geom_point(aes(color = as.factor(waifw_id))) + 
  labs(y = "S") +
  scale_color_brewer(palette = "Set1", labels = waifw_labs) +
  theme_bw()

jcm_all_v_long %>%
  filter(age < 100) %>%
  summarize(value = sum(value), .by = c("variable", "waifw_id", "v")) %>%
  filter(variable == "I") %>%
  ggplot(aes(x = v, y = value)) + 
  geom_line(aes(color = as.factor(waifw_id))) +
  geom_point(aes(color = as.factor(waifw_id))) + 
  theme_bw()
  
jcm_all_v_long %>%
  filter(variable == "S") %>%
  mutate(age = round(age, 3)) %>%
  left_join(data.frame(age = round(age_classes_jcm,3), 
                       bin_width = bin_width_jcm)) %>%
  ggplot(aes(x = age, y = value/bin_width, color = as.factor(v))) + 
  geom_line() +
  facet_wrap(vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) +
  scale_color_viridis_d() + 
  theme_bw()

jcm_all_v_long %>%
  filter(variable == "S") %>%
  mutate(age = round(age, 3)) %>%
  left_join(data.frame(age = round(age_classes_jcm,3), 
                       bin_width = bin_width_jcm)) %>%
  ggplot(aes(x = age, y = value/bin_width, color = as.factor(waifw_id))) + 
  geom_line() +
  facet_wrap(vars(v)) +
  scale_color_brewer(palette = "Set1", labels = waifw_labs) + 
  theme_bw()

jcm_all_v_long %>%
  filter(variable %in% c("S", "E", "I", "R")) %>%
  summarize(value = sum(value), .by = c("age", "v", "waifw_id")) %>%
  mutate(age = round(age, 3)) %>%
  left_join(data.frame(age = round(age_classes_jcm,3), 
                       bin_width = bin_width_jcm)) %>%
  filter(age < 100) %>%
  ggplot(aes(x = age, y = value/bin_width, color = as.factor(waifw_id))) + 
  geom_line() +
  facet_wrap(vars(v)) +
  scale_color_brewer(palette = "Set1", labels = waifw_labs) + 
  theme_bw()
    

start_vax = 0.94
release_vax = 0.6
chosen_dt = 1/52
# show equilibrium values for different vax rates and waifw matrices
jcm_release = vector("list", length(waifw))
# IC_release = vector("list", length(waifw))
for(i in 1:length(waifw)){
  print(paste0(i, "/", length(waifw)))
  # setup parameters
  paras_tmp = paras_jcm
  if(i > 1){
    s = chosen_scalars %>% filter(waifw_id == i, round(v,3) == start_vax) %>% pull(best_scalar)
    if(is.na(s)){print(paste0("NA SCALAR, SKIPPING (i = ", i)); next}
    paras_tmp["beta0"] = paras_tmp["beta0"] * s
  }
  IC_manual = jcm_all_v_long %>% filter(round(v,3) == start_vax, waifw_id == i, !(variable %in% c("BH", "BHs", "BHi", "BHb")))
  names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
  IC_manual = IC_manual$value
  names(IC_manual) = names_IC
  # try putting all I individuals in one compartment (say age 5)
  I_indx = which(substr(names(IC_manual), 1, 1) == "I")
  tot_I = sum(IC_manual[I_indx])
  IC_manual[I_indx] = 0
  IC_manual[which(names(IC_manual) == "I_5.5")] = 1
  IC_manual[which(names(IC_manual) == "R_5.5")] = IC_manual[which(names(IC_manual) == "R_10.5")] - (tot_I - 1)
  jcm_release[[i]] = run_ode(
    age_classes = age_classes_jcm, mort = mort,
    fert = fert, start_pop = paras_jcm["N"],
    compartments = compartments, dt = chosen_dt,
    vax_change_times = c(0), vax_rates = c(release_vax), waifw = waifw2[[i]],
    IC_type = "manual", IC_manual = IC_manual, max_t = 10, params = paras_tmp, #, IC_manual = new_IC
    adjust_beta_flag = FALSE, plot_flag = FALSE) %>%
    mutate(waifw_id = i, s = ifelse(i == 1, 1, s))
}
beep()

jcm_release_long = bind_rows(jcm_release)

unity_beta_long = jcm_release_long %>%
  mutate(s = ifelse(waifw_id == 1, 1, s)) %>% # for now to fix error
  filter(variable %in% c("C", "BHb", "BHs", "BHi")) %>%
  left_join(jcm_release_long %>%
              filter(variable %in% c("C", "BHb", "BHs", "BHi"), waifw_id == 1) %>%
              select(time, variable, value) %>%
              rename(homog_value = value)) %>%
  mutate(unity_beta = homog_value/(s*value))

Rt_long = jcm_release_long %>% 
  filter(variable == "S") %>%
  summarize(Rt = get_Rt(waifw2[[waifw_id]], value, paras_jcm["beta0"]*mean(s), paras_jcm["gamma"], paras_jcm["N"]), .by = c("time", "waifw_id"))

# try making plot with each waifw in one col (easier to interpret?)
pt0 = lapply(waifw2, melt) %>%
  bind_rows(.id = "waifw_id") %>%
  mutate(value_scaled = value/max(value), .by = c("waifw_id")) %>%
  ggplot(aes(x = Var1, y = Var2, fill = value_scaled)) + 
  geom_tile() + 
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) +
  scale_fill_distiller(palette = "YlGnBu") + 
  scale_x_continuous(expand = c(0,0),
                     name = "age (years)") +
  scale_y_continuous(expand = c(0,0),
                     name = "age (years)") +
  theme(legend.position = "none", 
        strip.background = element_blank())
pt1 = jcm_release_long %>%
  filter(variable %in% c("I")) %>%
  summarize(value = sum(value), .by = c("time", "variable", "waifw_id")) %>%
  ggplot(aes(x = time, y = value, color = as.factor(waifw_id))) + 
  geom_line(data = jcm_release_long %>%
              filter(variable %in% c("I"), waifw_id == 1) %>%
              summarize(value = sum(value), .by = c("time", "variable", "waifw_id")) %>% select(-waifw_id), linewidth = 0.8, color = "black") +
  geom_line(linewidth = 0.8) + 
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since release") + 
  scale_y_continuous(name = "infections") +
  theme_bw() + 
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor.y = element_blank(), 
        strip.background = element_blank())
pt2 = unity_beta_long %>%
  filter(variable == "C") %>%
  ggplot(aes(x = time, y = log(unity_beta), color = as.factor(waifw_id))) + 
  geom_text(data = data.frame(y = c(Inf, -Inf), x = c(10, 10), vjust = c(1, 0),
                              waifw_id = 1,
                              lab = c("\nslowing down\nwith age structure", "speeding up\nwith age structure\n")),
            aes(x = x, y = y, label = lab, vjust = vjust), hjust = 1, color = "black", size = 3, alpha = 1) +
  geom_line(data = unity_beta_long %>%
              filter(variable == "C", waifw_id == 1) %>% select(-waifw_id), linewidth = 0.8, color = "black") + 
  geom_line(aes(linetype = variable), linewidth = 0.8) + 
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  guides(color = FALSE) +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since release") + 
  scale_y_continuous(limits = c(-7.5, 7.5), name = "log(beta hat)") +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor.y = element_blank(), 
        strip.background = element_blank())
pt3 = Rt_long %>% 
  ggplot(aes(x = time, y = Rt, color = as.factor(waifw_id))) + 
  geom_hline(yintercept = 1, linetype = "dotted") + 
  geom_line(data = Rt_long %>% filter(waifw_id == 1) %>% select(-waifw_id), 
            linewidth = 0.8, color = "black") +
  geom_line(linewidth = 0.8) + 
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since release") + 
  scale_y_continuous(name = "Rt") +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor.y = element_blank(), 
        strip.background = element_blank())
pt4 = jcm_release_long %>% filter(variable %in% c("S", "I"), age <= 10) %>%
  mutate(bin_width = bin_width_jcm[as.factor(age)], 
         value = value/bin_width) %>% 
  dcast(time + age + waifw_id ~ variable, value.var = "value") %>%
  ggplot(aes(x = time, y = age)) + 
  geom_tile(aes(fill = S)) + 
  geom_point(data = jcm_release_long %>% filter(time %in% seq(0, 10, 0.1), variable == "I", age <= 10) %>%
               mutate(bin_width = bin_width_jcm[as.factor(age)],
                      I = ifelse(value < 1, NA, value/bin_width)),
             aes(size = I), color = "white") +
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) +
  scale_fill_viridis_c() + 
  scale_size_continuous(range = c(0,2)) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0,10,2), name = "years since release") + 
  scale_y_continuous(expand = c(0,0), breaks = seq(0,10,2), name = "age (year)") + 
  theme(legend.position = "none",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())
pt0/pt1/pt2/pt3/pt4

# show different beta hat values
unity_beta_long %>%
  # filter(variable == "C") %>%
  filter(waifw_id != 1) %>%
  ggplot(aes(x = time, y = log(unity_beta), color = as.factor(waifw_id))) + 
  geom_text(data = data.frame(y = c(Inf, -Inf), x = c(10, 10), vjust = c(1, 0),
                              waifw_id = 2,
                              lab = c("\nslowing down\nwith age structure", "speeding up\nwith age structure\n")),
            aes(x = x, y = y, label = lab, vjust = vjust), hjust = 1, color = "black", size = 3, alpha = 1) +
  geom_hline(yintercept = 0, linewidth = 0.8) +
  geom_line(aes(linetype = variable), linewidth = 0.8, alpha = 0.7) + 
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  guides(color = FALSE) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_linetype_manual(values = c("dotted", "twodash", "dashed", "solid"), labels = c("both homogeneous", "homogeneous I", "homogeneous S", "both structured")) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since release") + 
  scale_y_continuous(limits = c(-7.5, 7.5), name = "log(beta hat)") +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = "right",
        panel.grid.minor.y = element_blank(), 
        strip.background = element_blank())


pq0 = unity_beta_long %>%
  filter(variable == "C") %>%
  ggplot(aes(x = time, y = log(unity_beta), color = as.factor(waifw_id))) + 
  geom_line(data = unity_beta_long %>%
              filter(variable == "C", waifw_id == 1) %>% select(-waifw_id), linewidth = 0.8, color = "black") + 
  geom_line(linewidth = 0.8) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since release") + 
  scale_y_continuous(limits = c(-7.5, 7.5), name = "log(beta hat)") +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor.y = element_blank(), 
        strip.background = element_blank())
pq2 = jcm_release_long %>% filter(variable == "I") %>%
  mutate(age = round(age,3)) %>%
  left_join(data.frame(age = round(age_classes_jcm,3), 
                       bin_width = bin_width_jcm)) %>%
  # compute quantities of interest
  mutate(midpoint = age - bin_width/2) %>%
  mutate(mean_age = sum(value*midpoint)/sum(value), .by = c("time", "waifw_id")) %>%
  summarize(mean_age = mean(mean_age), 
            var_age = sum(value*(midpoint-mean_age)^2)/sum(value), 
            .by = c("time", "waifw_id")) %>%
  mutate(cv = sqrt(var_age)/mean_age) %>%
  melt(c("time" , "waifw_id")) %>%
  filter(variable != "var_age") %>%
  ggplot(aes(x = time, y = value, color = as.factor(waifw_id))) + 
  geom_line(size = 0.8) + 
  facet_grid(rows = vars(variable), scales = "free") +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since release") + 
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.title = element_blank())
pq0/pq2 + plot_layout(heights = c(0.33, 0.66))


#### RELEASE (WITH AGE STRUCTURE) FROM EVERY END VALUE -------------------------
chosen_dt = 1/52
chosen_start_vals = test_v[c(9)] #6, 8, 10
sim_all_release = vector("list", length(chosen_start_vals))
for(i in 1:length(chosen_start_vals)){
  print(paste0(i, "/", length(chosen_start_vals)))
  start_vax = chosen_start_vals[i]
  new_end_vals = test_v[c(1, 3, 5, 7)]
  new_end_vals = new_end_vals[new_end_vals < start_vax]
  tmp1 = vector("list", length(new_end_vals))
 for(j in 1:length(new_end_vals)){
   release_vax = new_end_vals[j]
   print(paste0(j, "/", length(new_end_vals)))
   tmp2 = vector("list", length(waifw2))
   for(k in 1:length(waifw2)){
     print(paste0(k, "/", length(waifw2)))
     paras_tmp = paras_jcm
     if(k > 1){
       s = chosen_scalars %>% filter(waifw_id == k, v == start_vax) %>% pull(best_scalar)
       if(is.na(s)){print(paste0("NA SCALAR, SKIPPING (k = ", k)); next}
       paras_tmp["beta0"] = paras_tmp["beta0"] * s
     }
     IC_manual = jcm_all_v_long %>% filter(v == start_vax, waifw_id == k, !(variable %in% c("BH", "BHs", "BHi", "BHb")))
     names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
     IC_manual = IC_manual$value
     names(IC_manual) = names_IC
     # try putting all I individuals in one compartment (say age 5)
     I_indx = which(substr(names(IC_manual), 1, 1) == "I")
     tot_I = sum(IC_manual[I_indx])
     IC_manual[I_indx] = 0
     IC_manual[which(names(IC_manual) == "I_5.5")] = 1
     IC_manual[which(names(IC_manual) == "R_5.5")] = IC_manual[which(names(IC_manual) == "R_5.5")] - (tot_I - 1)
     tmp2[[k]] = run_ode(
       age_classes = age_classes_jcm, mort = mort,
       fert = fert, start_pop = paras_jcm["N"],
       compartments = compartments, dt = chosen_dt,
       vax_change_times = c(0), vax_rates = c(release_vax), waifw = waifw2[[k]],
       IC_type = "manual", IC_manual = IC_manual, max_t = 8, params = paras_tmp, #, IC_manual = new_IC
       adjust_beta_flag = FALSE, plot_flag = FALSE) %>%
       mutate(waifw_id = k, s = ifelse(k == 1, 1, s), start_vax = start_vax, release_vax = release_vax)
   }
   tmp1[[j]] = bind_rows(tmp2)
 } 
  sim_all_release[[i]] = bind_rows(tmp1)
}
beep()

sim_all_release_long = bind_rows(sim_all_release)

unity_beta_long = sim_all_release_long %>%
  filter(variable %in% c("C", "BHb", "BHs", "BHi")) %>%
  left_join(sim_all_release_long %>%
              filter(variable %in% c("C", "BHb", "BHs", "BHi"), waifw_id == 1) %>%
              select(time, start_vax, release_vax, variable, value) %>%
              rename(homog_value = value)) %>%
  mutate(unity_beta = homog_value/(s*value))

Rt_long = sim_all_release_long %>% 
  filter(variable == "S") %>%
  summarize(Rt = get_Rt(waifw2[[waifw_id]], value, paras_jcm["beta0"]*mean(s), paras_jcm["gamma"], paras_jcm["N"]), .by = c("time", "waifw_id", "start_vax", "release_vax"))


p1 = sim_all_release_long %>%
  filter(variable == "I", release_vax %in% c(0, 0.3, 0.6, 0.9)) %>% 
  summarize(value = sum(value), .by = c("time", "waifw_id", "start_vax", "release_vax")) %>%
  ggplot(aes(x = time, y = value, color = as.factor(waifw_id), alpha = as.factor(-release_vax))) + 
  geom_line(linewidth = 0.8) + 
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  guides(color = FALSE) +
  scale_alpha_discrete(name = "release vax", range = c(0.4, 1)) +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw()
p2 = ggplot(data = unity_beta_long %>% filter(variable == "C", release_vax %in% c(0, 0.3, 0.6, 0.9)), 
       aes(x = time, y = log(unity_beta), color = as.factor(waifw_id), alpha = as.factor(-release_vax))) + 
  geom_hline(yintercept = 0, linewidth = 0.8) +
  geom_text(data = data.frame(y = c(Inf, -Inf), x = c(8,8), vjust = c(1, 0),
                              waifw_id = 1,
                              lab = c("\nslowing down\nwith age structure", "speeding up\nwith age structure\n")),
            aes(x = x, y = y, label = lab, vjust = vjust), hjust = 1, color = "black", size = 3, alpha = 1) +
  geom_line(linewidth = 0.8) + 
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  guides(color = FALSE) +
  scale_alpha_discrete(name = "release vax", range = c(0.4, 1)) +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_y_continuous(limits = c(-9, 9)) +
  theme_bw()
p3 = Rt_long %>% 
  filter(release_vax %in% c(0, 0.3, 0.6, 0.9)) %>%
  ggplot(aes(x = time, y = Rt, color = as.factor(waifw_id), alpha = as.factor(-release_vax))) + 
  geom_hline(yintercept = 1, linetype = "dotted") + 
  geom_line() + 
  guides(color = FALSE) +
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_alpha_discrete(name = "release vax", range = c(0.4, 1)) +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw()
p1/p2/p3

