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
source("R/twopatch-age-struc-functions.R")

#### MODEL SETUP ---------------------------------------------------------------
source("R/1_setup-WAIFW.R") # this will add waifw to environment, which is list of WAIFW matrices to test
n_waifw = length(waifw)

compartments = c("S1", "I1", "R1", "S2", "I2", "R2")

paras = c(mu = 1/50, N = 1, beta0 = 365, beta1 = 0, omega = 0,
              gamma = 365/14, pi = 1, delta = 0, ve = 0.95, p = 0)
paras["N"] = 500000
paras["delta"] = 1e-4

test_pi = c(0.2, 0.4, 0.6, 0.8, 1)

# R0
with(as.list(paras), beta0/(gamma + mu))

#### FIND EQUILIBRIUM ACROSS DIFFERENT VALUES OF PI ----------------------------
all_pi_steady = vector("list", length(test_pi))
for(j in 1:length(test_pi)){ #1
  print(paste0("j: ", j, "/", length(test_pi)))
  tmp = vector("list", length(waifw))
  for(i in 1:length(waifw)){
    print(paste0("i: ", i, "/", length(waifw)))
    IC = setup_IC(start_pop = paras_tmp["N"], age_classes, compartments, fert = fert, 
                  mort = mort, IC_type = "std") 
    Fmat <- buildFMatrix(age.classes = age_classes, fert = fert, ncompartments = length(compartments), time.step = 1)
    tmp_steady = runsteady(y = IC, times = c(0, 500), func = sirtwopatch_age_structured, parms = paras, # runsteady arguments
                           compartments = compartments, age_classes = age_classes, mort = mort, fert = fert,
                           vax_change_times = c(0), vax_rates = 0.94, waifw = waifw[[i]],
                           Fmat = Fmat, adjust_beta_flag = FALSE, print_warnings_flag = FALSE
    )
    tmp[[i]] = data.frame(variable = names(tmp_steady$y), 
                          value = tmp_steady$y, 
                          time = attributes(tmp_steady)$time,
                          pi = test_pi[j], waifw_id = i)
  }
  all_pi_steady[[j]] = bind_rows(tmp)
}
beep()

all_pi_steady_long = bind_rows(all_pi_steady) %>%
  mutate(
    age = ifelse(variable %in% c("C", "BH", "BHs", "BHi", "BHb"), NA, as.double(substr(variable, 4, nchar(as.character(variable))))),
    variable = ifelse(variable %in% c("C", "BH", "BHs", "BHi", "BHb"), variable, substr(variable, 1, 2)), 
    patch = ifelse(variable %in% c("C", "BH", "BHs", "BHi", "BHb"), NA, substr(variable, 2, 2))
  )
# add age- and patch-specific population sizes
all_pi_steady_long <- all_pi_steady_long %>%
  bind_rows(
    all_pi_steady_long %>% filter(!(variable %in% c("C", "BH", "BHs", "BHi", "BHb"))) %>% summarize(value = sum(value), .by = c("time", "waifw_id", "pi", "age", "patch")) %>% mutate(variable = paste0("N", patch))
  )

# and calculate R0
all_pi_steady_Rt_long = all_pi_steady_long %>% 
  filter(variable %in% c("S1", "S2", "N1", "N2")) %>%
  dcast(time + pi + waifw_id + age ~ variable) %>%
  summarize(Rt = get_Rt_twopatch(waifw = waifw[[waifw_id]], S1 = S1, S2 = S2, N1 = N1, N2 = N2, beta = paras["beta0"], gamma = paras["gamma"], pi = mean(pi)), 
                                 .by = c("time", "waifw_id", "pi"))

ggplot(data = all_pi_steady_Rt_long, aes(x = pi, y = Rt, color = as.factor(waifw_id))) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw()

ggplot(data = all_pi_steady_long %>% filter(nchar(variable) == 2) %>% mutate(variable = substr(variable, 1, 1)) %>%
         summarize(value = sum(value), .by = c("time", "pi", "waifw_id", "patch", "variable")), 
       aes(x = pi, y = value, color = as.factor(waifw_id))) + 
  geom_point() + 
  geom_line(aes(linetype = patch)) + 
  facet_wrap(vars(variable), scales = "free") + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw()

# for now, scalars to get to the same R0
scalars = all_pi_steady_Rt_long %>%
  left_join(all_pi_steady_Rt_long %>% filter(waifw_id == 1) %>% select(-waifw_id, -time) %>% rename(Rt_flat = Rt)) %>%
  mutate(scalar = Rt_flat/Rt)
  
#### FIND SCALARS TO YIELD EQUIVALENT EQUILIBRIUM INFECTIONS -------------------
# TO DO: test the R0 method here (will require fixing findStableAge)
#        find scalars across values of pi to test

#### SIMULATE TO VAX-ERA EQUILIBRIUM -------------------------------------------

#### LOOK AT REBOUND EFFECTS UNDER DIFFERENT VALUES OF PI ----------------------
chosen_dt = 1/52
# show equilibrium values for different vax rates and waifw matrices
twopatch_release = vector("list", length(waifw))
# IC_release = vector("list", length(waifw))
for(j in 1:length(test_pi)){
  tmp = vector("list", length(waifw))
  for(i in 1:length(waifw)){
    print(paste0(i, "/", length(waifw)))
    # setup parameters
    paras_tmp = paras
    paras_tmp["pi"] = test_pi[j]
    if(i > 1){
      s = scalars %>% filter(waifw_id == i, pi == test_pi[j]) %>% pull(scalar)
      if(is.na(s)){print(paste0("NA SCALAR, SKIPPING (i = ", i)); next}
      paras_tmp["beta0"] = paras_tmp["beta0"] * s
    }
    IC_manual = all_pi_steady_long %>% filter(waifw_id == i, pi == test_pi[j], !(variable %in% c("BH", "BHs", "BHi", "BHb", "N1", "N2", "C")))
    names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
    IC_manual = IC_manual$value
    names(IC_manual) = names_IC
    # try putting all I individuals in vax compartment (say age 5) first
    I1_indx = which(substr(names(IC_manual), 1, 1) == "I1")
    tot_I = sum(IC_manual[I1_indx])
    IC_manual[I1_indx] = 0
    IC_manual[which(names(IC_manual) == "I1_5.5")] = 1
    IC_manual[which(names(IC_manual) == "R1_5.5")] = IC_manual[which(names(IC_manual) == "R_5.5")] + (tot_I - 1)
    tmp[[i]] = run_ode(
      age_classes = age_classes, mort = mort,
      fert = fert, start_pop = paras["N"],
      compartments = compartments, dt = chosen_dt,
      vax_change_times = c(0), vax_rates = c(0.2), waifw = waifw[[i]],
      IC_type = "manual", IC_manual = IC_manual, max_t = 10, params = paras_tmp, #, IC_manual = new_IC
      adjust_beta_flag = FALSE, plot_flag = FALSE, func = sirtwopatch_age_structured) %>%
      mutate(waifw_id = i, s = ifelse(i == 1, 1, s), pi = test_pi[j])
  }
  twopatch_release[[j]] = bind_rows(tmp)
}
beep()

twopatch_release_long = bind_rows(twopatch_release) %>%
  mutate(
    patch = ifelse(variable %in% c("C", "BH", "BHs", "BHi", "BHb"), NA, substr(variable, 2, 2))
  )
# add age- and patch-specific population sizes
twopatch_release_long <- twopatch_release_long %>%
  bind_rows(
    twopatch_release_long %>% filter(!(variable %in% c("C", "BH", "BHs", "BHi", "BHb"))) %>% summarize(value = sum(value), .by = c("time", "waifw_id", "pi", "age", "patch")) %>% mutate(variable = paste0("N", patch))
  )

ggplot(data = twopatch_release_long %>% filter(substr(variable,1,1) == "I") %>% summarize(value = sum(value), .by = c('time', 'waifw_id', "pi")), 
       aes(x = time, y = value, color = as.factor(waifw_id), alpha = as.factor(pi))) + 
  geom_line() + 
  facet_wrap(vars(waifw_id)) +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw()

ggplot(data = twopatch_release_long %>% filter(substr(variable,1,1) == "I") %>% summarize(value = sum(value), .by = c('time', 'waifw_id', "pi", "patch")), 
       aes(x = time, y = value, color = as.factor(waifw_id), linetype = as.factor(patch))) + 
  geom_line() + 
  facet_grid(cols = vars(waifw_id), rows = vars(pi)) +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw()

twopatch_release_Rt_long = twopatch_release_long %>% 
  filter(variable %in% c("S1", "S2", "N1", "N2")) %>%
  dcast(time + pi + waifw_id + age ~ variable) %>%
  summarize(Rt = get_Rt_twopatch(waifw = waifw[[waifw_id]], S1 = S1, S2 = S2, N1 = N1, N2 = N2, beta = paras["beta0"], gamma = paras["gamma"], pi = mean(pi)), 
            .by = c("time", "waifw_id", "pi", ))

ggplot(data = twopatch_release_Rt_long, 
       aes(x = time, y = Rt, color = as.factor(waifw_id), alpha = as.factor(pi))) + 
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_line() + 
  facet_wrap(vars(waifw_id)) +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw()

# want to include beta hat? what does it mean here? probably could be extended...

#### RUN SIMULAITONS (OLD) -----------------------------------------------------
tst = vector("list", n_waifw)
for(i in 1:1){#n_waifw){
  print(paste0("i: ", i, "/", n_waifw))
  tst[[i]] = run_ode(
    age_classes = age_classes, mort = mort,
    fert = fert, start_pop = paras["N"],
    compartments = compartments, dt = 1/52,
    vax_change_times = c(0), vax_rates = c(0.5), waifw = waifw[[i]],
    IC_type = "std", max_t = 20, params = paras, func = sirtwopatch_age_structured,
    adjust_beta_flag = FALSE, plot_flag = FALSE)
}
beep()

tst_long = bind_rows(tst, .id = "waifw_id") %>%
  mutate(patch = ifelse(variable %in% c("C", "BH", "BHs", "BHi", "BHb"), NA, substr(variable, 2, 2)))
# add population size
tst_long = tst_long %>%
  bind_rows(tst_long %>% 
              filter(!(variable %in% c(c("C", "BH", "BHs", "BHi", "BHb")))) %>%
              summarize(value = sum(value), .by = c("waifw_id", "time", "age", "patch")) %>%
              mutate(variable = paste0("N", patch)))

p1 = ggplot(data = tst_long %>% 
         filter(substr(variable,1,1) == "I") %>% 
         summarize(value = sum(value), .by = c("time", "variable", "waifw_id")), 
       aes(x = time, y = value, color = variable)) + 
  geom_line() + 
  facet_wrap(vars(waifw_id))

ggplot(data = tst_long %>% 
         filter(substr(variable,1,1) == "N", variable != "NNA") %>% 
         summarize(value = sum(value), .by = c("time", "variable", "waifw_id")), 
       aes(x = time, y = value, color = variable)) + 
  geom_line() + 
  facet_wrap(vars(waifw_id))

p2 = tst_long %>% 
  filter(variable %in% c("S1", "S2", "N1", "N2")) %>%
  dcast(time + waifw_id + age ~ variable) %>%
  summarize(Rt = get_Rt_twopatch(waifw = waifw[[waifw_id]], S1 = S1, S2 = S2, N1 = N1, N2 = N2, beta = paras["beta0"], gamma = paras["gamma"], pi = paras["pi"]), 
            .by = c("time", "waifw_id")) %>%
  ggplot(aes(x = time, y = Rt)) + 
  geom_line()

p1/p2

