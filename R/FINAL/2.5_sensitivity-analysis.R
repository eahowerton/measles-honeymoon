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

tst_R0 = seq(13, 19, 4) #2
# tst_mu = 1/(seq(50, 80, 15))
tst_mu = seq(1/50, 1/80, length.out = 3) #6

paras_to_test = expand.grid(
  R0 = tst_R0, mu = tst_mu, waifw_id = 1:5, # parameters to change
  N = 500000, beta1 = 0, gamma = 365/14, p = 0, delta = 0)  %>% # other fixed parameters
  mutate(beta0 = (gamma + mu)*R0) # calculate beta to achieve same R0

start_vax = 0.95#seq(0.93, 0.97, 0.02)

#### GET R0 VALUES FOR EACH WAIFW ----------------------------------------------
stable_age = findStableStruct(age_classes, mort, fert, 1/52)$stable.age

find_scalar = function(s, R0, waifw, S, beta0, gamma, mu, N){
  # print(paste0("s: ", s, " R0: ", get_Rt(waifw, S, beta0*s, gamma, N)))
  diff = get_Rt(waifw, S, beta0*s, gamma, mu, N) - R0
  return(abs(diff))
}

paras_to_test$scalar_minimum = 0
paras_to_test$scalar_objective = 0

for(i in 1:nrow(paras_to_test)){
  o = optimize(f = find_scalar, tol = 1e-8, interval = c(0, 100), R0 = as.double(paras_to_test[i, "R0"]),
               waifw = waifw[[as.integer(paras_to_test[i, "waifw_id"])]], S = stable_age, 
               beta0 = as.double(paras_to_test[i, "beta0"]), gamma = as.double(paras_to_test[i, "gamma"]),  
               mu = as.double(paras_to_test[i, "mu"]), N = 1)
  paras_to_test[i, "scalar_minimum"] = o$minimum
  paras_to_test[i, "scalar_objective"] = o$objective
}

paras_to_test$beta0 = paras_to_test$beta0*paras_to_test$scalar_minimum

#### CONFIRM EQUILIBRIA ARE CLOSE, AFTER USING THESE SCALARS -------------------
vax_equilib = vector("list", length(waifw))
for(j in 1:length(start_vax)){
  print(paste0("j: ", j, "/", length(start_vax)))
  tmp = vector("list", length(waifw))
  for(i in 1:nrow(paras_to_test)){
    print(i)
    Fmat = buildFMatrix(age.classes = age_classes, fert = fert, ncompartments = length(compartments))
    IC = setup_IC(start_pop = unlist(paras_to_test[i, "N"]), age_classes = age_classes, 
                  compartments = compartments, mort = mort, fert = fert, 
                  IC_type = "stable-age-noI", dt = 1/52)
    t = runsteady(
      y = IC, times = c(0, 750), func = sir_age_structured, parms = unlist(paras_to_test[i,]), # runsteady arguments
      compartments = compartments, age_classes = age_classes, mort = mort, fert = fert,
      vax_change_times = c(0), vax_rates = c(start_vax[j]),
      waifw = waifw[[as.integer(paras_to_test[i, "waifw_id"])]], length(age_classes), rtol = 1e-12, atol = 1e-12,
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


#### CALCULATE HONEYMOON TIME --------------------------------------------------
# i.e., simulate across the full range of release vax values
release_vax_full = seq(0, 0.9, 0.1)
nyears_postrelease = 20
chosen_dt = 1/12

# generate data.frame of simulations to run
# note: do not need to repeat for all WAIFWs because there are no new infections
tst_release_full = expand.grid(start_vax = start_vax, 
                               release_vax = release_vax_full, 
                               para_id = 1:nrow(paras_to_test %>% select(R0, mu, N, beta1, gamma, p, delta) %>% unique())) %>%
  left_join(paras_to_test %>% select(R0, mu, N, beta1, gamma, p, delta) %>% unique() %>% mutate(para_id = seq_len(n()))) %>%
  mutate(beta0 = R0*(gamma + mu))

susc_after_release_full = vector("list", nrow(tst_release_full))
for(i in 1:nrow(tst_release_full)){
  print(paste0("i: ", i, "/", nrow(tst_release_full)))
  tmp_start_vax = unlist(tst_release_full[i, "start_vax"])
  tmp_release_vax = unlist(tst_release_full[i, "release_vax"])
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
    age_classes = age_classes, mort = mort, fert = fert, start_pop = unlist(tst_release_full[i, "N"]),
    compartments = compartments, dt = chosen_dt, waifw = waifw[[1]],
    vax_change_times = c(0), vax_rates = c(tmp_release_vax), 
    IC_type = "manual", IC_manual = IC_manual, max_t = nyears_postrelease, 
    params = unlist(tst_release_full[i, ]),
    adjust_beta_flag = FALSE, plot_flag = FALSE
  )
}
beep()

susc_after_release_long_full = bind_rows(susc_after_release_full, .id = "sim_id") %>%
  mutate(sim_id = as.integer(sim_id)) %>%
  left_join(tst_release_full %>% mutate(sim_id = seq_len(dplyr::n())))

rt_after_release_full = expand.grid(para_id = 1:nrow(tst_release_full), 
                                    waifw_id = 1:5) %>% 
  left_join(tst_release_full %>% mutate(para_id = seq_len(n()))) %>%
  select(-beta0) %>% 
  left_join(paras_to_test)

rt_after_release_rslt = vector("list", nrow(rt_after_release_full))
for(i in 1:nrow(rt_after_release_full)){
  print(paste0("i: ", i, "/", nrow(rt_after_release_full)))
  paras_tmp = unlist(rt_after_release_full[i, ])
  rt_after_release_rslt[[i]] = susc_after_release_long_full %>% 
    filter(variable == "S", R0 == unlist(rt_after_release_full[i, "R0"]),
           mu == unlist(rt_after_release_full[i, "mu"]), 
           start_vax == unlist(rt_after_release_full[i, "start_vax"]), 
           release_vax == unlist(rt_after_release_full[i, "release_vax"]) ) %>%
    summarize(Rt = get_Rt(waifw = waifw[[as.integer(rt_after_release_full[i, "waifw_id"])]], 
                          S = value, beta0 = paras_tmp["beta0"], gamma = paras_tmp["gamma"], 
                          N = paras_tmp["N"], mu = paras_tmp["mu"]), 
              .by = c("time", "start_vax", "release_vax", "R0", "mu"))
}
beepr::beep()

rt_after_release_full_long = bind_rows(rt_after_release_rslt, .id = "par_id") %>%
  left_join(rt_after_release_full %>% mutate(par_id = seq_len(n())))

saveRDS(rt_after_release_full_long, "R/FINAL/data/rt_after_release_sensitivity_long.rds")

honeymoon_period = rt_after_release_full_long %>%
  filter(Rt > 1) %>%
  mutate(min_time = min(time), .by = c("waifw_id", "start_vax", "release_vax", "R0", "mu")) %>%
  filter(time == min_time) %>%
  select(-min_time)

honeymoon_period %>%
  ggplot(aes(x = release_vax, y = time, color = as.factor(waifw_id))) +
  geom_line(linewidth = 0.8) + 
  facet_grid(rows = vars(mu), cols = vars(R0, start_vax)) +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0, 9, 0.3), label = scales::percent, name = "coverage after drop") + 
  scale_y_continuous(name = "time to Re > 1") + 
  theme_bw() + 
  theme(legend.title = element_blank())

ggplot(honeymoon_period %>% filter(release_vax > 0.6, R0 == 17, round(start_vax,3) == 0.95, waifw_id == 5), aes(x = mu, y = release_vax, fill = as.factor(waifw_id), alpha = time)) + 
  geom_tile() + 
  geom_contour(aes(z = time), color = "black") +
  facet_grid(rows = vars(R0, start_vax), cols = vars(waifw_id)) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank())
  

ggplot(honeymoon_period %>% filter(release_vax > 0.6, R0 == 17, round(start_vax,3) == 0.95, waifw_id == 5)) + 
  geom_line(aes(x = mu, y = time)) + 
  facet_wrap(vars(release_vax), scales = "free")

