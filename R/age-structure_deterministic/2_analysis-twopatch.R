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

source("R/age-structure_deterministic/0_helper-functions.R")
source("R/age-structure_deterministic/0_SIR-age.R")
source("R/twopatch-age-struc-functions.R")

#### MODEL SETUP ---------------------------------------------------------------
source("R/age-structure_deterministic/1_setup-WAIFW.R") # this will add waifw to environment, which is list of WAIFW matrices to test
n_waifw = length(waifw)

compartments = c("S", "I", "R")
compartments_twopatch = c("S1", "I1", "R1", "S2", "I2", "R2")

R0 <- 14
paras = c(mu = 1/50, N = 500000, beta1 = 0,
          gamma = 365/14, delta = 1e-4, p = 0, phi = 1) # delta = 1e-4
paras["beta0"]  = (paras["gamma"] + paras["mu"])*R0

paras_noimport = paras
paras_noimport["delta"] = 0

# R0
with(as.list(paras), beta0/(gamma + mu))

prop_N1 = 0.9
N1 = prop_N1*paras["N"]
N2 = (1-prop_N1)*paras["N"]
# fix population 1 vaccination rates, and vary population 2 vax rates for different overall coverage
start_vax = 0.96#c(0.9, 0.92, 0.94, 0.96)
start_vax1 = rep(0.96, length(start_vax))
start_vax2 = (start_vax*paras["N"] - N1*start_vax1)/N2 
release_vax1 = seq(0.7, 0.9, 0.02)
release_vax2 = seq(0, 0.8, 0.02)

tst_phi = seq(0, 1, 0.5)

stable_age = findStableStruct(age_classes, mort, fert, 1/52)$stable.age

#### SHOW HOW R0 DEPENDS ON V1, V2, AND N1/(N1 + N2) ---------------------------
tst_all_paras = expand.grid(prop_N1 = seq(0.1,0.9,0.1), 
                            v1 = seq(0, 1, 0.1),
                            v2 = seq(0, 1, 0.1),
                            waifw_id = 1:5, 
                            phi = seq(0, 1, 0.1))

tst_all_paras = tst_all_paras %>%
  mutate(row_id = seq_len(dplyr::n())) %>%
  mutate(R0 = get_Rt_twopatch(waifw = waifw[[waifw_id]], 
                              S1 = stable_age*prop_N1*(1-v1), 
                              S2 = stable_age*(1-prop_N1)*(1-v2), 
                              N1 = prop_N1, 
                              N2 = (1-prop_N1), 
                              beta0 = paras["beta0"], gamma = paras["gamma"], mu = paras["mu"], phi = phi), 
         .by = c("row_id"))

ggplot(data = tst_all_paras %>% mutate(max_R0 = max(R0), .by = c("waifw_id")) %>% filter(waifw_id == 3), 
       aes(x = v1, y = v2)) + 
  geom_tile(aes(fill = R0)) + 
  geom_contour(aes(z = R0)) +
  facet_grid(rows = vars(prop_N1), cols = vars(phi)) + 
  scale_fill_viridis_c()

ggplot(data = tst_all_paras %>% 
         filter(waifw_id == 3, phi %in% seq(0, 1, 0.5), round(prop_N1,4) %in% c(0.2, 0.5, 0.8)), 
       aes(x = v1, y = v2)) + 
  geom_tile(aes(fill = R0)) + 
  geom_contour(aes(z = R0), breaks = 10, color = "black") +
  facet_grid(rows = vars(prop_N1), cols = vars(phi), labeller = label_both) + 
  scale_fill_viridis_c()


# get contours
contour_df = tst_all_paras %>%
  filter(phi %in% seq(0, 1, 0.5), round(prop_N1,4) %in% c(0.2, 0.5, 0.8)) %>%
  arrange(v1, v2) 

contours = vector("list")
cntr = 1
for(i in unique(contour_df$waifw_id)){
  for(j in unique(contour_df$phi)){
    for(k in unique(contour_df$prop_N1)){
      print(cntr)
      tmp = contour_df %>% filter(waifw_id == i, phi == j, prop_N1 == k)
      tmp_x = tmp$v1 %>% unique()
      tmp_y = tmp$v2 %>% unique()
      tmp_z = contour_df %>% 
        filter(waifw_id == i, phi == j, prop_N1 == k) %>%
        dcast(v1 ~ v2, value.var = "R0") %>%
        select(-v1) %>%
        as.matrix()
      if(max(tmp_z)< 10){
        contours[[cntr]] = data.frame(
          waifw_id = i, phi = j, prop_N1 = k, 
          v1 = tmp2[[1]]$x, 
          v2 = tmp2[[1]]$y,
          R0 = NA
        )
        break
      }
      tmp2 = contourLines(tmp_x, tmp_y, tmp_z, levels = c(10))
      contours[[cntr]] = data.frame(
        waifw_id = i, phi = j, prop_N1 = k, 
        v1 = tmp2[[1]]$x, 
        v2 = tmp2[[1]]$y,
        R0 = tmp2[[1]]$level
      )
      cntr = cntr + 1
    }
  }
}

contours = bind_rows(contours)

ggplot(data = contours, 
       aes(x = v1, y = v2, color = phi, linetype = as.factor(prop_N1), group = interaction(prop_N1, phi))) + 
  geom_path() +
  facet_wrap(vars(waifw_id))

# now try analytical result
tst_again = tst_all_paras %>% 
  filter(prop_N1 == 0.5) %>%
  mutate(P1_waifw0 = get_Rt(waifw[[waifw_id]], 
                            S = stable_age*(1-v1)*(prop_N1), 
                            beta0 = 1, gamma = 1, mu = 0, N = (prop_N1)), 
         P2_waifw0 = get_Rt(waifw[[waifw_id]], 
                            S = stable_age*(1-v2), 
                            beta0 = 1, gamma = 1, mu = 0, N = 1),
        P1_R0 = get_Rt(waifw[[waifw_id]], 
                        S = stable_age*(1-v1), 
                        beta0 = paras["beta0"], gamma = paras["gamma"], mu = paras["mu"], N = 1), 
         P2_R0 = get_Rt(waifw[[waifw_id]], 
                        S = stable_age*(1-v2), 
                        beta0 = paras["beta0"], gamma = paras["gamma"], mu = paras["mu"], N = 1), 
         .by = "row_id") %>%
  left_join(data.frame(waifw_id = 1:5) %>%
              mutate(waifw0 = get_Rt(waifw[[waifw_id]], S = stable_age, beta0 = 1, gamma = 1, mu = 0, N = 1), .by = c("waifw_id")))

tst_again = tst_again %>%
  mutate(R0_analyt =  paras["beta0"]/( paras["gamma"] +  paras["mu"]) * (P1_waifw0*P2_waifw0*(1-phi^2))/(phi*P1_waifw0 - P2_waifw0)) %>%
  mutate(R0_analyt2 = ((1-v1)*(1-v2)*(1-phi^2))/(phi*(1-v1) - (1-v2))*waifw0*paras["beta0"]/(paras["gamma"] + paras["mu"]))

ggplot(data = tst_again %>% filter(v2 == 0.5, phi %in% c(0,0.5,1), prop_N1 == 0.5), 
       aes(x = v1)) + 
  geom_line(aes(y = R0, color = "sim")) + 
  geom_line(aes(y = R0_analyt, color = "analyt1")) + 
  geom_line(aes(y = R0_analyt2, color = "analyt12")) + 
  facet_grid(cols = vars(waifw_id), rows = vars(phi)) +
  theme_bw()
  

# let's dive in for waifw_id == 2, v1 == 0.6, v2 == 0.5, phi == 0
tst_again %>% filter(waifw_id == 1, phi == 0, v2 == 0.5)

get_Rt_twopatch(waifw[[1]], S1 = stable_age, S2 = stable_age, v1 == 0.8, v2 == 0.5,  beta0 = paras["beta0"], gamma = paras["gamma"], mu = paras["mu"], phi == 0)

#### GET R0 VALUES FOR EACH WAIFW ----------------------------------------------

find_scalar = function(s, R0, waifw, S1, S2, N1, N2, beta0, gamma, mu, phi){
  diff = get_Rt_twopatch(waifw, S1, S2, N1, N2, beta0*s, gamma, mu, phi) - R0
  return(abs(diff))
}

scalars = expand.grid(waifw_id = 1:5, phi = tst_phi, scalar = NA, diff = NA)
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
paras_all_noimport = lapply(1:nrow(scalars), function(i){paras_tmp = paras_noimport; paras_tmp["beta0"] = paras_tmp["beta0"]*scalars[i,3]; return(paras_tmp)})


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
                                   release_vax2 = release_vax2) %>%
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
    params = paras_all_noimport[[1]], # can use any because differences don't matter without infections (here only tracking S)
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
rt_after_release_twopatch = vector("list", nrow(scalars))
p_extinct_after_release_twopatch = vector("list", nrow(scalars))
for(i in 1:nrow(scalars)){
  print(paste0(i, "/", nrow(scalars)))
  paras_tmp = paras_all_noimport[[i]]
  tmp_waifw_id = scalars[i, "waifw_id"]
  tmp_phi = scalars[i, "phi"]
  rt_after_release_twopatch[[i]] = susc_after_release_twopatch_wide %>% 
    filter(time < 10) %>%
    summarize(Rt = get_Rt_twopatch(waifw = waifw[[tmp_waifw_id]], S1 = S1, S2 = S2, N1 = N1, N2 = N2, 
                                   beta0 = paras_tmp["beta0"], gamma = paras_tmp["gamma"], mu = paras_tmp["mu"], phi = tmp_phi), 
              .by = c("time", "start_vax", "release_vax1", "release_vax2"))
  # p_extinct_after_release_twopatch[[i]] = susc_after_release_twopatch_wide %>% 
  #   filter(time < 10) %>%
  #   dplyr::reframe(compute_extinction_prob_twopatch(
  #     waifw = waifw[[tmp_waifw_id]], S1 = S1, S2 = S2, N1 = N1, N2 = N2, 
  #     beta0 = paras_tmp["beta0"], gamma = paras_tmp["gamma"], mu = paras_tmp["mu"], phi = tmp_phi, age_classes),
  #     .by = c("time", "start_vax", "release_vax1", "release_vax2")) 
}

rt_after_release_twopatch_long = bind_rows(rt_after_release_twopatch, .id = "row_id") %>%
  mutate(row_id = as.integer(row_id)) %>%
  left_join(scalars %>% select(waifw_id, phi) %>% mutate(row_id = seq_len(n())))
saveRDS(rt_after_release_twopatch_long, "data/rt_after_release_twopatch_long.rds")
# p_extinct_after_release_twopatch_long = bind_rows(p_extinct_after_release_twopatch, .id = "row_id") %>%
#   mutate(row_id = as.integer(row_id)) %>%
#   left_join(scalars %>% select(waifw_id, phi) %>% mutate(row_id = seq_len(n())))

# some plots to explore the effects of phi
ggplot(data = rt_after_release_twopatch_long %>% filter(round(release_vax1,3) == 0.9, release_vax2 == 0),
       aes(x = time, y = Rt, color = as.factor(phi))) +
  geom_hline(yintercept = 1) +
  geom_line() +
  facet_wrap(vars(waifw_id), labeller = labeller(waifw_id = waifw_labs), scales = "free") +
  theme_bw()

susc_after_release_twopatch_wide %>% 
  filter(time < 10, round(release_vax1,3) == 0.9, release_vax2 == 0) %>%
  summarize(Rt = get_Rt_twopatch(waifw = waifw[[2]], S1 = S1, S2 = S2, N1 = N1, N2 = N2, 
                                 beta0 = paras_tmp["beta0"], gamma = paras_tmp["gamma"], mu = paras_tmp["mu"], phi = tmp_phi), 
            .by = c("time", "start_vax", "release_vax1", "release_vax2")) %>%
  ggplot(aes(x = time, y = Rt, color = phi))
  

# honeymmon
honeymoon_period = rt_after_release_twopatch_long %>%
  filter(Rt > 1.05) %>%
  mutate(min_time = min(time), .by = c("waifw_id", "phi", "start_vax", "release_vax1", "release_vax2")) %>%
  filter(time == min_time) %>%
  select(-min_time)

honeymoon_period %>% filter(start_vax == 0.96) %>%
  mutate(release_coverage = (release_vax1*N1 + release_vax2*N2)/paras["N"]) %>%
  ggplot(aes(x = release_vax1, y = release_vax2)) + 
  geom_tile(aes(fill = time)) + 
  geom_contour(aes(z = release_coverage), color = "gray") + 
  geom_contour(aes(z = time), color = "red") +
  # metR::geom_text_contour(aes(label = time), color = "white") +
  facet_grid(cols = vars(waifw_id), rows = vars(phi)) + 
  scale_fill_viridis_c()

honeymoon_period %>% filter(start_vax == 0.96) %>%
  ggplot(aes(x = release_vax1, y = release_vax2)) + 
  geom_contour(aes(z = time, color = as.factor(waifw_id)), breaks = c(2, 5, 7)) + 
  facet_grid(cols = vars(phi)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.position = "bottom")

honeymoon_period %>% filter(start_vax == 0.96) %>%
  ggplot(aes(x = phi, y = time, color = as.factor(waifw_id))) +
  geom_point() +
  geom_line(linewidth = 0.8) +
  facet_grid(cols = vars(release_vax1), rows = vars(release_vax2)) +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.position = "bottom")


honeymoon_period %>% filter(start_vax == 0.96, round(release_vax1,3) %in% c(0.7, 0.8, 0.9)) %>%
  ggplot(aes(x = release_vax2, y = time, color = as.factor(waifw_id), alpha = as.factor(phi))) +
  geom_line(linewidth = 0.8) +
  facet_grid(cols = vars(release_vax1)) +
  scale_alpha_manual(values = c(1, 0.7, 0.4)) +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.position = "bottom")

p_extinct_after_release_twopatch_long %>%
  filter(age < 15, time < 12, start_vax == 0.96) %>%
  ggplot(aes(x = time, y = age))+ 
  geom_tile(aes(fill = outbreak_prob)) + 
  geom_contour(aes(z = outbreak_prob), color = "white") +
  facet_grid(rows = vars(waifw_id, patch), cols = vars(phi), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_fill_viridis_c() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0),
                     # breaks = which(age_classes %in% c(2, 4, 6, 8, 10, 30, 50, 70)),
                     # labels = c(2, 4, 6, 8, 10, 30, 50, 70),
                     name = "age (years)") +
  theme_bw()

#### try simulating some releases ----------------------------------------------
some_twopatch_releases = vector("list", nrow(tst_release_twopatch))
for(j in 1:length(waifw)){
  print(paste0(j, "/", length(waifw)))
  tmp = vector("list", length(tst_phi))
  for(i in 1:length(tst_phi)){
    tmp_paras = paras_all_noimport[[j]]
    tmp_paras["phi"] = tst_phi[i]
    tmp_start_vax = unlist(tst_release_twopatch[1, "start_vax"])
    tmp_release_vax1 = unlist(tst_release_twopatch[1, "release_vax1"])
    tmp_release_vax2 = unlist(tst_release_twopatch[1, "release_vax2"])
    IC_manual = vax_equilib_twopatch_long %>% filter(start_vax == tmp_start_vax, waifw_id == j, !(variable %in% c("C", "BH", "BHs", "BHi", "BHb")))
    names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
    IC_manual = IC_manual$value
    names(IC_manual) = names_IC
    # now move all infections into R (local extinction)
    I_indx = which(substr(names(IC_manual), 1, 1) == "I")
    R_indx = which(substr(names(IC_manual), 1, 1) == "R")
    IC_manual[R_indx] = IC_manual[R_indx] + IC_manual[I_indx]
    IC_manual[I_indx] = 0
    IC_manual["I1_5"] = IC_manual["I1_5"] + 1
    IC_manual["R1_5"] = IC_manual["R1_5"] - 1
    tmp[[i]] = run_ode_twopatch(
      age_classes = age_classes, mort = mort, fert = fert, start_pop = paras["N"],
      compartments = compartments_twopatch, dt = chosen_dt, waifw = waifw[[j]],
      vax_change_times = c(0), vax_rates1 = c(tmp_release_vax1), vax_rates2 = c(tmp_release_vax2),  
      IC_type = "manual", IC_manual = IC_manual, max_t = 15, 
      params = tmp_paras,
      adjust_beta_flag = FALSE, plot_flag = FALSE
    ) %>%
      mutate(phi = tst_phi[i])
  }
  some_twopatch_releases[[j]] = bind_rows(tmp)
}
beep()

some_twopatch_releases_long = bind_rows(some_twopatch_releases, .id = "waifw_id") %>%
  mutate(patch = substr(variable, 2,2), 
         compartment = substr(variable, 1, 1))

ggplot(data = some_twopatch_releases_long %>% filter(compartment == "I", time < 5) %>% summarize(value = sum(value), .by = c("time", "phi", "waifw_id", "patch")), 
       aes(x = time, y = value, color = as.factor(phi), linetype = as.factor(patch))) + 
  geom_line() + 
  facet_wrap(vars(waifw_id), scales = "free", labeller = labeller(waifw_id = waifw_labs)) + 
  theme_bw()

