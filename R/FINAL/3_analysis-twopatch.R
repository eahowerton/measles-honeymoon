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
source("R/FINAL/0_twopatch-age-functions.R")

#### MODEL SETUP ---------------------------------------------------------------
source("R/FINAL/2.1_setup-WAIFW.R") # this will add waifw to environment, which is list of WAIFW matrices to test
n_waifw = length(waifw)

compartments = c("S", "I", "R")
compartments_twopatch = c("S1", "I1", "R1", "S2", "I2", "R2")

R0 <- 17
paras = c(mu = 1/80, N = 500000, beta1 = 0,
          gamma = 365/14, delta = 0, p = 0, phi = 1) # delta = 1e-4
paras["beta0"]  = (paras["gamma"] + paras["mu"])*R0

start_vax = 0.95

prop_N1 = 0.9
N1 = prop_N1*paras["N"]
N2 = (1-prop_N1)*paras["N"]
# fix population 1 vaccination rates, and vary population 2 vax rates for different overall coverage
start_vax = 0.96
start_vax1 = start_vax
start_vax2 = (start_vax*paras["N"] - N1*start_vax1)/N2 
release_vax1 = 0.9
release_vax2 = seq(0, 0.8, 0.02)

tst_phi = seq(0, 1, 0.5)

stable_age = findStableStruct(age_classes, mort, fert, 1/52)$stable.age

#### GET R0 VALUES FOR EACH WAIFW ----------------------------------------------

find_scalar = function(s, R0, waifw, S1, S2, N1, N2, beta0, gamma, mu, phi){
  diff = get_Rt_twopatch(waifw, S1, S2, N1, N2, beta0*s, gamma, mu, phi) - R0
  return(abs(diff))
}

# fix R0 when phi = 1 (close, but not equivalent to single patch model)
scalars = expand.grid(waifw_id = 1:5, phi = 1, scalar = NA, diff = NA)
for(i in 1:nrow(scalars)){
  tmp_waifw_id = unlist(scalars[i, "waifw_id"])
  tmp_phi = unlist(scalars[i, "phi"])
  o = optimize(f = find_scalar, tol = 1e-8, interval = c(0, 100), R0 = R0, waifw = waifw[[tmp_waifw_id]],
               S1 = stable_age*prop_N1, S2 = stable_age*(1-prop_N1),
               N1 = prop_N1, N2 = (1-prop_N1),
               beta0 = paras["beta0"], gamma = paras["gamma"], mu = paras["mu"], phi = tmp_phi)
  scalars[i, "scalar"] = o$minimum
  scalars[i, "diff"] = o$objective
}

# create a list of parameters for 
paras_all = lapply(1:nrow(scalars), function(i){paras_tmp = paras; paras_tmp["beta0"] = paras_tmp["beta0"]*scalars[i,3]; return(paras_tmp)})

#### FIND PRE-RELEASE EQUILIBRIUM ----------------------------------------------
# equilibrium pre-release population size
vax_equilib_twopatch = vector("list", length(waifw))
for(j in 1:length(start_vax)){
  print(paste0("j: ", j, "/", length(start_vax)))
  tmp = vector("list", length(waifw))
  for(i in 1:length(waifw)){
    print(i)
    Fmat = buildFMatrix_twopatch(age.classes = age_classes, fert = fert, 
                                 ncompartments = length(compartments_twopatch), time.step = 1)
    IC = setup_IC(start_pop = paras["N"], age_classes = age_classes, 
                  compartments = compartments_twopatch, mort = mort, fert = fert, 
                  IC_type = "stable-age", dt = 1/52, rel_size = c(prop_N1, 1-prop_N1))
    t = runsteady(
      y = IC, times = c(0, 500), stol = 1e-10, func = sirtwopatch_age_structured, parms = paras_all[[i]], # run steady arguments
      compartments = compartments_twopatch, age_classes = age_classes, mort = mort, fert = fert,
      vax_change_times = c(0), vax_rates1 = c(start_vax1[j]), vax_rates2 = c(start_vax2[j]),
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
  vax_equilib_twopatch[[j]] = bind_rows(tmp)
}

vax_equilib_twopatch_long = bind_rows(vax_equilib_twopatch) %>%
  mutate(age = ifelse(substr(variable, 1, 1) %in% compartments, 
                      as.double(substr(variable, 4, nchar(variable))), NA), 
         patch = ifelse(substr(variable, 1, 1) %in% compartments, 
                        as.integer(substr(variable, 2, 2)), NA),
         variable = ifelse(substr(variable, 1, 1) %in% compartments, substr(variable, 1, 2), variable),
         variable_short = ifelse(substr(variable, 1, 1) %in% compartments, substr(variable, 1, 1), variable)) %>%
  mutate(age = round(age, 4)) %>%
  left_join(data.frame(age = round(age_classes, 4), 
                       bin_width = bin_width))

# note: do not need to repeat for all WAIFWs or values of phi
# because there are no new infections

#### RELEASE UNDER DIFFERENT COVERAGE VALUES FOR EACH PATCH
tst_release_twopatch = expand_grid(start_vax = start_vax, 
                                   release_vax1 = release_vax1, 
                                   release_vax2 = release_vax2
                                   ) %>%
  mutate(release_vax = (release_vax1*N1 + release_vax2*N2)/paras["N"])
chosen_dt = 1/26

susc_after_release_twopatch = vector("list", nrow(tst_release_twopatch))
for(i in 1:nrow(tst_release_twopatch)){
  print(paste0("i: ", i, "/", nrow(tst_release_twopatch)))
  tmp_start_vax = unlist(tst_release_twopatch[i, "start_vax"])
  tmp_release_vax1 = unlist(tst_release_twopatch[i, "release_vax1"])
  tmp_release_vax2 = unlist(tst_release_twopatch[i, "release_vax2"])
  IC_manual = vax_equilib_twopatch_long %>% filter(start_vax == tmp_start_vax, waifw_id == 1, !(variable %in% c("C", "BH", "BHs", "BHi", "BHb")))
  names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
  IC_manual = IC_manual$value
  names(IC_manual) = names_IC
  # now move all infections into R (local extinction)
  I_indx = which(substr(names(IC_manual), 1, 1) == "I")
  R_indx = which(substr(names(IC_manual), 1, 1) == "R")
  IC_manual[R_indx] = IC_manual[R_indx] + IC_manual[I_indx]
  IC_manual[I_indx] = 0
  susc_after_release_twopatch[[i]] = run_ode_twopatch(
    age_classes = age_classes, mort = mort, fert = fert, start_pop = paras["N"],
    compartments = compartments_twopatch, dt = chosen_dt, waifw = waifw[[1]],
    vax_change_times = c(0), vax_rates1 = c(tmp_release_vax1), vax_rates2 = c(tmp_release_vax2),  
    IC_type = "manual", IC_manual = IC_manual, max_t = 15, 
    params = paras_all[[1]], # can use any because differences don't matter without infections (here only tracking S)
    adjust_beta_flag = FALSE, plot_flag = FALSE
  )
}
beep()

susc_after_release_twopatch_long = bind_rows(susc_after_release_twopatch, .id = "sim_id") %>%
  mutate(sim_id = as.integer(sim_id)) %>%
  left_join(tst_release_twopatch %>% mutate(sim_id = seq_len(dplyr::n())))

susc_after_release_twopatch_wide = susc_after_release_twopatch_long %>%
  filter(substr(variable, 1, 1) == "S") %>%
  left_join(susc_after_release_twopatch_long %>% 
              filter(substr(variable, 1, 1) %in% c("S", "I", "R")) %>%
              mutate(patch = substr(variable, 2,2)) %>%
              summarize(value = sum(value), .by = c("sim_id", "time", "patch", "start_vax", "release_vax1", "release_vax2")) %>%
              mutate(variable = paste0("N", patch)) %>% 
              select(-patch) %>%
              dcast(sim_id + time + start_vax + release_vax1 + release_vax2 ~ variable, value.var = "value")) %>%
  dcast(sim_id + time + age + start_vax + release_vax1 + release_vax2 + N1 + N2 ~ variable, value.var = "value")

# now repeat for each waifw
rt_release = expand.grid(waifw_id = 1:length(waifw), 
                         phi = c(0, 0.5, 1))
rt_after_release_twopatch = vector("list", nrow(rt_release))
p_extinct_after_release_twopatch = vector("list", nrow(rt_release))
for(i in 1:nrow(rt_release)){
  print(paste0(i, "/", nrow(rt_release)))
  tmp_waifw_id = rt_release[i, "waifw_id"]
  tmp_phi = rt_release[i, "phi"]
  paras_tmp = paras_all[[tmp_waifw_id]]
  rt_after_release_twopatch[[i]] = susc_after_release_twopatch_wide %>% 
    filter(time < 10) %>%
    summarize(Rt = get_Rt_twopatch(waifw = waifw[[tmp_waifw_id]], S1 = S1, S2 = S2, N1 = N1, N2 = N2, 
                                   beta0 = paras_tmp["beta0"], gamma = paras_tmp["gamma"], mu = paras_tmp["mu"], phi = tmp_phi), 
              .by = c("time", "start_vax", "release_vax1", "release_vax2"))
  p_extinct_after_release_twopatch[[i]] = susc_after_release_twopatch_wide %>%
    filter(time < 10) %>%
    dplyr::reframe(compute_extinction_prob_twopatch(
      waifw = waifw[[tmp_waifw_id]], S1 = S1, S2 = S2, N1 = N1, N2 = N2,
      beta0 = paras_tmp["beta0"], gamma = paras_tmp["gamma"], mu = paras_tmp["mu"], phi = tmp_phi, age_classes),
      .by = c("time", "start_vax", "release_vax1", "release_vax2"))
}

rt_after_release_twopatch_long = bind_rows(rt_after_release_twopatch, .id = "row_id") %>%
  mutate(row_id = as.integer(row_id)) %>%
  left_join(rt_release %>% select(waifw_id, phi) %>% mutate(row_id = seq_len(n())))
# saveRDS(rt_after_release_twopatch_long, "R/FINAL/data/rt_after_release_twopatch_long.rds") # rt_after_release_twopatch_long_fixphi1.rds is the final version (for now)
# rt_after_release_twopatch_long_diffscalar = readRDS("data/rt_after_release_twopatch_long_fixphi1.rds")


# honeymmon
honeymoon_period = rt_after_release_twopatch_long %>%
  filter(Rt > 1) %>%
  summarize(min_time = min(time), .by = c("waifw_id", "phi", "start_vax", "release_vax1", "release_vax2"))

honeymoon_period %>%
  ggplot(aes(x = release_vax2, y = min_time, color = as.factor(waifw_id), alpha = as.factor(phi))) +
  geom_line(linewidth = 0.8) +
  scale_alpha_manual(values = rev(c(1, 0.7, 0.4)), name = "patch mixing") +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs, name = "") +
  scale_x_continuous(labels = scales::percent, name = "release vax % (patch 2)") +
  scale_y_continuous(name = "time to Re > 1 (years)") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("R/FINAL/figures/honeymoon_age_twopatch.pdf", width = 6, height = 5)

# probability of outbreak by age and time
p_extinct_after_release_twopatch_long = bind_rows(p_extinct_after_release_twopatch, .id = "row_id") %>%
  mutate(row_id = as.integer(row_id)) %>%
  left_join(rt_release %>% select(waifw_id, phi) %>% mutate(row_id = seq_len(n())))
# saveRDS(p_extinct_after_release_twopatch_long, "R/FINAL/data/p_extinct_after_release_twopatch_long.rds") # rt_after_release_twopatch_long_fixphi1.rds is the final version (for now)

p_extinct_after_release_twopatch_long %>% filter(release_vax2 == 0.6, age < 15) %>%
  mutate(patch = paste0("P", patch)) %>%
  dcast(start_vax + release_vax1 + release_vax2 + age + waifw_id + phi + time ~ patch, value.var = "outbreak_prob") %>%
  ggplot(aes(x = time, y = age)) + 
  geom_tile(aes(fill = P1-P2)) + 
  facet_grid(cols = vars(waifw_id), rows = vars(phi), labeller = label_both) + 
  scale_fill_viridis_c(option = "mako") + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_bw()
ggsave("R/FINAL/figures/p_outbreak.pdf", width = 8, height = 4)

