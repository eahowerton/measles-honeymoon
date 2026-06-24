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
library(stringr)
library(readxl)

source("R/0_helper-functions.R")
source("R/0_SIR-age-functions.R")

#### MODEL SETUP ---------------------------------------------------------------
source("R/2.1_setup-WAIFW.R") # this will add waifw to environment, which is list of WAIFW matrices to test
n_waifw = length(waifw)

compartments = c("S", "I", "R")

R0 <- 17
paras = c(mu = 1/80, N = 500000, beta1 = 0,
          gamma = 365/14, delta = 0, p = 0, phi = 1) # delta = 1e-4
paras["beta0"]  = (paras["gamma"] + paras["mu"])*R0

tst_mu = seq(-5.5, -3.5, length.out = 20)#seq(log(1/45), log(1/200), length.out = 20)

start_vax = 0.95

#### GET R0 VALUES FOR EACH WAIFW ----------------------------------------------
scalars_all = vector("list", length(tst_mu))
vax_equilib_long = vector("list", length(tst_mu))
for(j in 1:length(tst_mu)){
  tmp_mu = exp(tst_mu[j])
  fert = rep(tmp_mu, length(age_classes))
  mort = rep(tmp_mu, length(age_classes))
  stable_age = findStableStruct(age_classes, mort, fert, 1/52)$stable.age
  scalars_tmp = data.frame(waifw_id = 5, scalar = NA, diff = NA)
  i = 5
  # for(i in 5:length(waifw)){
  o = optimize(f = find_scalar, tol = 1e-8, interval = c(0, 100), R0 = R0,
               waifw = waifw[[5]], S = stable_age, beta0 = paras["beta0"], gamma = paras["gamma"],  mu = tmp_mu, N = 1)
  print(get_Rt(waifw[[5]], stable_age, paras["beta0"]*o$minimum, paras["gamma"],  mu = tmp_mu, N = 1))
  scalars_tmp[1, 2:3] = c(o$minimum, o$objective)
  # scalars_tmp[i, 2:3] = c(o$minimum, o$objective)
  # }
  scalars_all[[j]]  = scalars_tmp
  # pre-vax equilibrium = (1-vax_cov)*N(a) where N(a) is stable age distribution
  vax_equilib_long[[j]] = expand.grid(age = age_classes, 
                                      start_vax = start_vax, 
                                      variable = c("S", "I", "R")) %>%
    left_join(data.frame(age = age_classes, 
                         bin_width = bin_width, 
                         N = stable_age*paras["N"])) %>%
    mutate(value = ifelse(variable == "S", N*(1-start_vax), ifelse(variable == "I", 0, N*start_vax)))
}

# create a list of parameters for 
paras_all = lapply(1:length(tst_mu), function(i){paras_tmp = paras; paras_tmp["beta0"] = paras_tmp["beta0"]*scalars_all[[i]][,2]; paras_tmp["mu"] = exp(tst_mu[i]); return(paras_tmp)})

#### CALCULATE HONEYMOON TIME --------------------------------------------------
# i.e., simulate across the full range of release vax values
release_vax_full = seq(0.5, 1, 0.02)
nyears_postrelease = 20
chosen_dt = 1/52

# generate data.frame of simulations to run
# note: do not need to repeat for all WAIFWs because there are no new infections
susc_after_release_full = vector("list", length(tst_mu))
for(j in 1:length(tst_mu)){
  print(paste0("j: ", j, "/", length(tst_mu)))
  tmp_mu = exp(tst_mu[j])
  fert = rep(tmp_mu, length(age_classes))
  mort = rep(tmp_mu, length(age_classes))
  tst_release_full = expand.grid(start_vax = start_vax, 
                                 release_vax = release_vax_full)
  tmp = vector("list", nrow(tst_release_full))
  for(i in 1:nrow(tst_release_full)){
    tmp_start_vax = tst_release_full[i, "start_vax"]
    tmp_release_vax = tst_release_full[i, "release_vax"]
    IC_manual = vax_equilib_long[[j]] %>% filter(start_vax == tmp_start_vax) %>% arrange(age)
    names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
    IC_manual = IC_manual$value
    names(IC_manual) = names_IC
    tmp[[i]] = run_ode(
      age_classes = age_classes, mort = mort, fert = fert, start_pop = paras["N"],
      compartments = compartments, dt = chosen_dt, waifw = waifw[[5]],
      vax_rates = c(tmp_release_vax), 
      IC_manual = IC_manual, max_t = nyears_postrelease, 
      params = paras_all[[j]],
      adjust_beta_flag = FALSE
    )
  }
  susc_after_release_full[[j]] = bind_rows(tmp, .id = "sim_id") %>%
    mutate(sim_id = as.integer(sim_id)) %>%
    left_join(tst_release_full %>% mutate(sim_id = seq_len(dplyr::n())))
}
beep()

susc_after_release_long_full = bind_rows(susc_after_release_full, .id = "mu_id") %>%
  mutate(mu_id = as.integer(mu_id)) %>%
  left_join(bind_rows(lapply(paras_all, function(i){data.frame(t(i))}) ) %>% mutate(mu_id = seq_len(length(tst_mu)))) 

# POLYMOD ONLY
rt_after_release_full = susc_after_release_long_full %>% 
  filter(variable == "S") %>%
  summarize(Rt = get_Rt(waifw = waifw[[5]], S = value, beta0 = beta0, gamma = gamma, N = N, mu = mu), 
            .by = c("time", "start_vax", "release_vax", "mu_id", "mu"))

honeymoon_period = rt_after_release_full %>%
  left_join(bind_rows(lapply(paras_all, function(i){data.frame(t(i))}) ) %>% mutate(mu_id = seq_len(length(tst_mu)))) %>%
  filter(Rt > 1) %>%
  mutate(min_time = min(time), .by = c("mu_id", "start_vax", "release_vax", "mu")) %>%
  filter(time == min_time) %>%
  select(-min_time)

saveRDS(honeymoon_period, "data/output-data/honeymoon_period_by_birthrate.rds")

