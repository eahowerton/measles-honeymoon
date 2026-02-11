library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidytable)
# library(epimdr2)
library(cowplot)
library(dplyr)
library(beepr)
library(fields)
library(rootSolve)
library(patchwork)

source("R/FINAL/0_helper-functions.R")
source("R/FINAL/0_SIR-age-functions.R")

#### MODEL SETUP ---------------------------------------------------------------
source("R/FINAL/2.1_setup-WAIFW.R") # this will add waifw to environment, which is list of WAIFW matrices to test
n_waifw = length(waifw)

compartments = c("S", "I", "R")

R0 <- 17
paras = c(mu = 1/80, N = 500000, beta1 = 0,
          gamma = 365/14, delta = 0, p = 0, phi = 1) # delta = 1e-4
paras["beta0"]  = (paras["gamma"] + paras["mu"])*R0

fert = rep(paras["mu"], length(age_classes))
mort = rep(paras["mu"], length(age_classes))
# mort <- c(rep(0,length(age_classes)-1),1)
# fert <-  c(rep(0,length(age_classes)-1),1)

start_vax = 0.96

#### GET R0 VALUES FOR EACH WAIFW ----------------------------------------------
stable_age = findStableStruct(age_classes, mort, fert, 1/52)$stable.age

find_scalar = function(s, R0, waifw, S, beta0, gamma, mu, N){
  # print(paste0("s: ", s, " R0: ", get_Rt(waifw, S, beta0*s, gamma, N)))
  diff = get_Rt(waifw, S, beta0*s, gamma, mu, N) - R0
  return(abs(diff))
}

scalars = data.frame(waifw_id = 1:5, scalar = NA, diff = NA)
for(i in 1:length(waifw)){
  o = optimize(f = find_scalar, tol = 1e-8, interval = c(0, 100), R0 = R0,
               waifw = waifw[[i]], S = stable_age, beta0 = paras["beta0"], gamma = paras["gamma"],  mu = paras["mu"], N = 1)
  print(get_Rt(waifw[[i]], stable_age, paras["beta0"]*o$minimum, paras["gamma"],  mu = paras["mu"], N = 1))
  scalars[i, 2:3] = c(o$minimum, o$objective)
}

# create a list of parameters for 
paras_all = lapply(1:5, function(i){paras_tmp = paras; paras_tmp["beta0"] = paras_tmp["beta0"]*scalars[i,2]; return(paras_tmp)})
paras_all_noimport = lapply(1:5, function(i){paras_tmp = paras; paras_tmp["beta0"] = paras_tmp["beta0"]*scalars[i,2]; return(paras_tmp)})

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
                  IC_type = "stable-age-noI", dt = 1/52)
    t = runsteady(
      y = IC, times = c(0, 750), func = sir_age_structured, parms = paras_all[[i]], # runsteady arguments
      compartments = compartments, age_classes = age_classes, mort = mort, fert = fert,
      vax_change_times = c(0), vax_rates = c(start_vax[j]),
      waifw = waifw[[i]], length(age_classes), rtol = 1e-12, atol = 1e-12,
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
  filter(variable == "S") %>%
  mutate(value_scaled = value/bin_width) %>%
  ggplot(aes(x = age, y = value_scaled, color = as.factor(waifw_id))) +
  geom_line() + 
  geom_line(data = data.frame(x = age_classes, y = stable_age*paras["N"]*(1-start_vax)/(bin_width)),
  aes(x = x, y = y), linetype = "dotted", color = "black") +
  facet_wrap(vars(start_vax), scales = 'free') +
  theme_bw()

vax_equilib_long %>%
  tidyr::spread(variable, value) %>%
  mutate(N = S + I + R, prop_S = S/N, prop_age = N/paras["N"]) %>%
  filter(waifw_id == 1) %>%
  ggplot(aes(x = age, y = prop_age/bin_width)) + 
  geom_line() + 
  geom_line(data = data.frame(x = age_classes, y = stable_age/bin_width),
            aes(x = x, y = y), linetype = "dotted", color = "black")

vax_equilib_long %>%
  

#### CALCULATE HONEYMOON TIME --------------------------------------------------
# i.e., simulate across the full range of release vax values
release_vax_full = seq(0, 0.9, 0.01)
nyears_postrelease = 20
chosen_dt = 1/52

# generate data.frame of simulations to run
# note: do not need to repeat for all WAIFWs because there are no new infections
tst_release_full = expand.grid(start_vax = start_vax, 
                               release_vax = release_vax_full)

susc_after_release_full = vector("list", nrow(tst_release_full))
for(i in 1:nrow(tst_release_full)){
  print(paste0("i: ", i, "/", nrow(tst_release_full)))
  tmp_start_vax = tst_release_full[i, "start_vax"]
  tmp_release_vax = tst_release_full[i, "release_vax"]
  IC_manual = vax_equilib_long %>% filter(start_vax == tmp_start_vax, waifw_id == 1, !(variable %in% c("C", "BH", "BHs", "BHi", "BHb")))
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
    summarize(Rt = get_Rt(waifw = waifw[[i]], S = value, beta0 = paras_tmp["beta0"], gamma = paras_tmp["gamma"], N = paras_tmp["N"], mu = paras_tmp["mu"]), 
              .by = c("time", "start_vax", "release_vax"))
}

rt_after_release_full_long = bind_rows(rt_after_release_full, .id = "waifw_id")

# save results
saveRDS(susc_after_release_long_full, "R/FINAL/data/susc_after_release_long_start96.rds")
saveRDS(rt_after_release_full_long, "R/FINAL/data/rt_after_release_long_start96.rds")

honeymoon_period = rt_after_release_full_long %>%
  filter(Rt > 1) %>%
  mutate(min_time = min(time), .by = c("waifw_id", "start_vax", "release_vax")) %>%
  filter(time == min_time) %>%
  select(-min_time)

honeymoon_period %>%
  ggplot(aes(x = release_vax, y = time, color = as.factor(waifw_id))) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0, 9, 0.3), label = scales::percent, name = "coverage after drop") +
  scale_y_continuous(name = "time to Re > 1") +
  theme_bw() +
  theme(legend.title = element_blank())
# ggsave("R/FINAL/figures/honeymoon_age_full_start95.pdf", width = 5, height = 4)


honeymoon_period96 = honeymoon_period


susc_after_release_long_full95 = readRDS("R/FINAL/data/susc_after_release_long.rds")
rt_after_release_full_long95 = readRDS("R/FINAL/data/rt_after_release_long.rds")


honeymoon_period95 = rt_after_release_full_long95 %>%
  filter(Rt > 1) %>%
  mutate(min_time = min(time), .by = c("waifw_id", "start_vax", "release_vax")) %>%
  filter(time == min_time) %>%
  select(-min_time)

honeymoon_period = honeymoon_period96 %>%
  bind_rows(honeymoon_period95)

honeymoon_period %>%
  ggplot(aes(x = release_vax, y = time, color = as.factor(waifw_id), alpha = as.factor(start_vax))) +
  geom_line(linewidth = 0.8) +
  facet_wrap(vars(waifw_id), labeller = labeller(waifw_id = waifw_labs), nrow = 1) + 
  guides(color = "none") + 
  scale_alpha_manual(name = "coverage before drop", values = c(0.6, 1), labels = c("95%", "96%")) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0, 9, 0.3), label = scales::percent, name = "coverage after drop") +
  scale_y_continuous(name = "time to Re > 1") +
  theme_bw() +
  theme(strip.background = element_blank(), 
        legend.position = "bottom", 
        panel.grid.minor = element_blank())
ggsave("R/FINAL/figures/honeymoon_age_full_comparestart.pdf", width = 9, height = 3)

