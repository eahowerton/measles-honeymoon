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
library(rootSolve)

source("R/age-structure_deterministic/0_helper-functions.R")
source("R/age-structure_deterministic/0_SIR-age.R")

#### MODEL SETUP ---------------------------------------------------------------
source("R/age-structure_deterministic/1_setup-WAIFW.R") # this will add waifw to environment, which is list of WAIFW matrices to test
n_waifw = length(waifw)

compartments = c("S", "I", "R")

R0 <- 14
paras = c(mu = 1/50, N = 500000, beta1 = 0,
          gamma = 365/14, delta = 1e-4, p = 0) # delta = 1e-4
paras["beta0"]  = (paras["gamma"] + paras["mu"])*R0

paras_noimport = paras
paras_noimport["delta"] = 0

# R0
with(as.list(paras), beta0/(gamma + mu))

start_vax = c(0.9, 0.92, 0.94, 0.96)

#### GET R0 VALUES FOR EACH WAIFW ----------------------------------------------
stable_age = findStableStruct(age_classes, mort, fert, 1/52)$stable.age

find_scalar = function(s, R0, waifw, S, beta0, gamma, N){
  # print(paste0("s: ", s, " R0: ", get_Rt(waifw, S, beta0*s, gamma, N)))
  diff = get_Rt(waifw, S, beta0*s, gamma, N) - R0
  return(abs(diff))
}

scalars = data.frame(waifw_id = 1:5, scalar = NA, diff = NA)
for(i in 1:length(waifw)){
  o = optimize(f = find_scalar, tol = 1e-8, interval = c(0, 100), R0 = R0,
               waifw = waifw[[i]], S = stable_age, beta0 = paras["beta0"], gamma = paras["gamma"], N = 1)
  print(get_Rt(waifw[[i]], stable_age, paras["beta0"]*o$minimum, paras["gamma"], 1))
  scalars[i, 2:3] = c(o$minimum, o$objective)
}

# create a list of parameters for 
paras_all = lapply(1:5, function(i){paras_tmp = paras; paras_tmp["beta0"] = paras_tmp["beta0"]*scalars[i,2]; return(paras_tmp)})
paras_all_noimport = lapply(1:5, function(i){paras_tmp = paras_noimport; paras_tmp["beta0"] = paras_tmp["beta0"]*scalars[i,2]; return(paras_tmp)})

#### CONFIRM EQUILIBRIA ARE CLOSE, AFTER USING THESE SCALARS -------------------
vax_equilib = vector("list", length(waifw))
for(j in 1:length(start_vax)){
  print(paste0("j: ", j, "/", length(start_vax)))
  tmp = vector("list", length(waifw))
  for(i in 1:length(waifw)){
    print(i)
    Fmat = buildFMatrix(age.classes = age_classes, fert = fert, ncompartments = length(compartments))
    IC = setup_IC(start_pop = paras["N"], age_classes = age_classes, 
                  compartments = compartments, mort = mort, fert = fert, 
                  IC_type = "stable-age", dt = 1/52)
    t = runsteady(
      y = IC, times = c(0, 750), func = sir_age_structured, parms = paras_all[[i]], # runsteady arguments
      compartments = compartments, age_classes = age_classes, mort = mort, fert = fert,
      vax_change_times = c(0), vax_rates = c(start_vax[j]),
      waifw = waifw[[i]], length(age_classes), 
      Fmat = Fmat, adjust_beta_flag = FALSE, print_warnings_flag = FALSE
      )
    tmp[[i]] = data.frame(variable = names(t$y), 
                          value = t$y, 
                          time = attr(t, "time"),
                          steady = attr(t, "steady"), 
                          waifw_id = i, 
                          start_vax = start_vax[j])
  }
  vax_equilib[[j]] = bind_rows(tmp)
}

vax_equilib_long = bind_rows(vax_equilib) %>%
  mutate(age = ifelse(substr(variable, 1, 1) %in% compartments, 
                      as.double(substr(variable, 3, nchar(variable))), NA), 
         variable = ifelse(substr(variable, 1, 1) %in% compartments, substr(variable, 1, 1), variable)) %>%
  mutate(age = round(age, 4)) %>%
  left_join(data.frame(age = round(age_classes, 4), 
                       bin_width = bin_width))

# plot results to do some checks
vax_equilib_long %>%
  filter(variable %in% c("S", "I", "R")) %>%
  summarize(value = sum(value), .by = c("waifw_id", "start_vax", "variable")) %>%
  ggplot(aes(x = start_vax, y = value, color = as.factor(waifw_id))) + 
  geom_point() + 
  geom_line() +
  facet_wrap(vars(variable), scales = "free") + 
  theme_bw()
  
vax_equilib_long %>%
  filter(variable == "I") %>%
  mutate(value_scaled = value/bin_width) %>%
  ggplot(aes(x = age, y = value_scaled, color = as.factor(waifw_id))) +
  geom_line() + 
  facet_wrap(vars(start_vax), scales = 'free') +
  theme_bw()

#### SIMULATE VACCINATION RELEASE FROM EQUILIBRIUM -----------------------------
# 1. remove all infections from the population (assume local extinction)
# 2. simulate accumulation of susceptibles over time at different release vax levels
release_vax = c(0, 0.4, 0.8)
nyears_postrelease = 12
chosen_dt = 1/52

# generate data.frame of simulations to run
# note: do not need to repeat for all WAIFWs because there are no new infections
tst_release = expand.grid(start_vax = start_vax, 
                         release_vax = release_vax)

susc_after_release = vector("list", nrow(tst_release))
for(i in 1:nrow(tst_release)){
  print(paste0("i: ", i, "/", nrow(tst_release)))
  tmp_start_vax = tst_release[i, "start_vax"]
  tmp_release_vax = tst_release[i, "release_vax"]
  IC_manual = vax_equilib_long %>% filter(start_vax == tmp_start_vax, waifw_id == 1, !(variable %in% c("BH", "BHs", "BHi", "BHb")))
  names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
  IC_manual = IC_manual$value
  names(IC_manual) = names_IC
  # now move all infections into R (local extinction)
  I_indx = which(substr(names(IC_manual), 1, 1) == "I")
  R_indx = which(substr(names(IC_manual), 1, 1) == "R")
  IC_manual[R_indx] = IC_manual[R_indx] + IC_manual[I_indx]
  IC_manual[I_indx] = 0
  susc_after_release[[i]] = run_ode(
    age_classes = age_classes, mort = mort, fert = fert, start_pop = paras["N"],
    compartments = compartments, dt = chosen_dt, waifw = waifw[[1]],
    vax_change_times = c(0), vax_rates = c(tmp_release_vax), 
    IC_type = "manual", IC_manual = IC_manual, max_t = nyears_postrelease, 
    params = paras_all_noimport[[1]],
    adjust_beta_flag = FALSE, plot_flag = FALSE
  )
}
beep()

susc_after_release_long = bind_rows(susc_after_release, .id = "sim_id") %>%
  mutate(sim_id = as.integer(sim_id)) %>%
  left_join(tst_release %>% mutate(sim_id = seq_len(dplyr::n())))

susc_after_release_long %>% 
  filter(variable == "S") %>%
  summarize(value = sum(value), .by = c("time", "start_vax", "release_vax")) %>%
  ggplot(aes(x = time, y = value, color = as.factor(release_vax))) + 
  geom_line() + 
  facet_wrap(vars(start_vax)) + 
  theme_bw()


#### CALCULATE RT AND PROBABILITY OF OUTBREAK ----------------------------------
# now repeat for each waifw
rt_after_release = vector("list", length(waifw))
p_extinct_after_release = vector("list", length(waifw))
for(i in 1:length(waifw)){
  print(i)
  paras_tmp = paras_all_noimport[[i]]
  rt_after_release[[i]] = susc_after_release_long %>% 
    filter(variable == "S") %>%
    summarize(Rt = get_Rt(waifw[[i]], value, paras_tmp["beta0"], paras_tmp["gamma"], paras_tmp["N"]), 
              .by = c("time", "start_vax", "release_vax"))
  p_extinct_after_release[[i]] = susc_after_release_long %>% 
    filter(variable == "S") %>%
    dplyr::reframe(compute_extinction_prob(waifw[[i]], value, paras_tmp["beta0"], paras_tmp["gamma"], paras_tmp["N"], age_classes), 
              .by = c("time", "start_vax", "release_vax"))
}

# let's plot and see what it looks like
rt_after_release_long = bind_rows(rt_after_release, .id = "waifw_id")
p_extinct_after_release_long = bind_rows(p_extinct_after_release, .id = "waifw_id")

ggplot(data = rt_after_release_long, 
       aes(x = time, y = Rt, color = as.factor(waifw_id))) + 
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_line() + 
  facet_grid(cols = vars(start_vax), rows = vars(release_vax), labeller = label_both) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.position = "bottom")

ggplot(data = p_extinct_after_release_long %>% filter(age == 5), 
       aes(x = time, y = outbreak_prob, color = as.factor(waifw_id))) + 
  geom_line() + 
  facet_grid(cols = vars(start_vax), rows = vars(release_vax), labeller = label_both) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.position = "bottom")

# choose highest contact group (for something more comparable?)
highest_contact = unlist(lapply(waifw, function(i){which.max(colSums(i))}))
highest_contact[1] = which(age_classes == 5) # modify flat to age 5 instead of infant
highest_contact = age_classes[highest_contact]

ggplot(data = p_extinct_after_release_long %>%
         mutate(highest_contact_flag = ifelse(age == highest_contact[waifw_id], TRUE, FALSE)) %>%
         filter(highest_contact_flag == TRUE),
       aes(x = time, y = outbreak_prob, color = as.factor(waifw_id))) + 
  geom_line() + 
  facet_grid(cols = vars(start_vax), rows = vars(release_vax), labeller = label_both) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.position = "bottom")

# put them together
p1 = ggplot(data = rt_after_release_long %>% filter(start_vax == 0.96, release_vax == 0.4), 
            aes(x = time, y = Rt, color = as.factor(waifw_id))) + 
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_line() + 
  facet_grid(cols = vars(start_vax), rows = vars(release_vax), labeller = label_both) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.position = "bottom")
p2 = ggplot(data = p_extinct_after_release_long %>% filter(age == 7, start_vax == 0.96, release_vax == 0.4), 
       aes(x = time, y = outbreak_prob, color = as.factor(waifw_id))) + 
  geom_line() + 
  facet_grid(cols = vars(start_vax), rows = vars(release_vax), labeller = label_both) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.position = "bottom")
p1/p2

#### FULL HONEYMOON TIME FIGURE ------------------------------------------------
# i.e., simulate across the full range of release vax values
release_vax_full = seq(0, 0.88, 0.04)
nyears_postrelease = 20
chosen_dt = 1/52

# generate data.frame of simulations to run
# note: do not need to repeat for all WAIFWs because there are no new infections
tst_release_full = expand.grid(start_vax = 0.94, 
                          release_vax = release_vax_full)

susc_after_release_full = vector("list", nrow(tst_release_full))
for(i in 1:nrow(tst_release_full)){
  print(paste0("i: ", i, "/", nrow(tst_release_full)))
  tmp_start_vax = tst_release_full[i, "start_vax"]
  tmp_release_vax = tst_release_full[i, "release_vax"]
  IC_manual = vax_equilib_long %>% filter(start_vax == tmp_start_vax, waifw_id == 1, !(variable %in% c("BH", "BHs", "BHi", "BHb")))
  names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
  IC_manual = IC_manual$value
  names(IC_manual) = names_IC
  # now move all infections into R (local extinction)
  I_indx = which(substr(names(IC_manual), 1, 1) == "I")
  R_indx = which(substr(names(IC_manual), 1, 1) == "R")
  IC_manual[R_indx] = IC_manual[R_indx] + IC_manual[I_indx]
  IC_manual[I_indx] = 0
  susc_after_release_full[[i]] = run_ode(
    age_classes = age_classes, mort = mort, fert = fert, start_pop = paras["N"],
    compartments = compartments, dt = chosen_dt, waifw = waifw[[1]],
    vax_change_times = c(0), vax_rates = c(tmp_release_vax), 
    IC_type = "manual", IC_manual = IC_manual, max_t = nyears_postrelease, 
    params = paras_all_noimport[[1]],
    adjust_beta_flag = FALSE, plot_flag = FALSE
  )
}
beep()

susc_after_release_long_full = bind_rows(susc_after_release_full, .id = "sim_id") %>%
  mutate(sim_id = as.integer(sim_id)) %>%
  left_join(tst_release_full %>% mutate(sim_id = seq_len(dplyr::n())))

rt_after_release_full = vector("list", length(waifw))
p_extinct_after_release_full = vector("list", length(waifw))
for(i in 1:length(waifw)){
  print(i)
  paras_tmp = paras_all_noimport[[i]]
  rt_after_release_full[[i]] = susc_after_release_long_full %>% 
    filter(variable == "S") %>%
    summarize(Rt = get_Rt(waifw[[i]], value, paras_tmp["beta0"], paras_tmp["gamma"], paras_tmp["N"]), 
              .by = c("time", "start_vax", "release_vax"))
  p_extinct_after_release_full[[i]] = susc_after_release_long_full %>% 
    filter(variable == "S") %>%
    dplyr::reframe(compute_extinction_prob(waifw[[i]], value, paras_tmp["beta0"], paras_tmp["gamma"], paras_tmp["N"], age_classes), 
                   .by = c("time", "start_vax", "release_vax"))
}

rt_after_release_full_long = bind_rows(rt_after_release_full, .id = "waifw_id")
p_extinct_after_release_full_long = bind_rows(p_extinct_after_release_full, .id = "waifw_id")

honeymoon_period = rt_after_release_full_long %>%
  filter(Rt > 1) %>%
  mutate(min_time = min(time), .by = c("waifw_id", "start_vax", "release_vax")) %>%
  filter(time == min_time) %>%
  select(-min_time)

honeymoon_period %>%
  ggplot(aes(x = release_vax, y = time, color = as.factor(waifw_id))) +
  geom_line(linewidth = 0.8) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.position = "bottom")

p_extinct_after_release_full_long %>%
  mutate(highest_contact_flag = ifelse(age == highest_contact[waifw_id], TRUE, FALSE)) %>%
  filter(highest_contact_flag == TRUE, release_vax %in% seq(0, 0.84, 0.12)) %>%
  ggplot(aes(x = time, y = outbreak_prob, color = as.factor(waifw_id), alpha = release_vax, group = interaction(waifw_id, release_vax))) + 
  geom_line() + 
  facet_wrap(vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  guides(color = FALSE) +
  scale_alpha_continuous(range = c(1, 0.3)) +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())

# show full heatmap of age-specific prob of outbreak 6 months after honeymoon period
p_extinct_after_release_full_long %>%
  left_join(honeymoon_period %>% rename(honeymoon_time = time)) %>%
  filter(release_vax %in% seq(0, 0.84, 0.12), time > honeymoon_time + 0.5) %>%
  mutate(min_time = min(time), .by = c("waifw_id", "start_vax", "release_vax")) %>%
  filter(time == min_time) %>%
  ggplot(aes(x = release_vax, y = as.factor(age), fill = outbreak_prob))+ 
  geom_tile() + 
  facet_wrap(vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_fill_viridis_c() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0),
                     # breaks = which(age_classes %in% c(2, 4, 6, 8, 10, 30, 50, 70)),
                     # labels = c(2, 4, 6, 8, 10, 30, 50, 70),
                     name = "age (years)") +
  theme_bw()

p_extinct_after_release_full_long %>%
  left_join(honeymoon_period %>% rename(honeymoon_time = time)) %>%
  filter(release_vax == 0.6, age < 15, time < 12) %>%
  ggplot(aes(x = time, y = age))+ 
  geom_tile(aes(fill = outbreak_prob)) + 
  geom_contour(aes(z = outbreak_prob), color = "white") +
  # geom_vline(aes(xintercept = honeymoon_time), color = "white") +
  facet_wrap(vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_fill_viridis_c() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0),
                   # breaks = which(age_classes %in% c(2, 4, 6, 8, 10, 30, 50, 70)),
                   # labels = c(2, 4, 6, 8, 10, 30, 50, 70),
                   name = "age (years)") +
  theme_bw()

#### TWO PATCH -----------------------------------------------------------------

