library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidytable)
library(cowplot)
library(dplyr)
library(beepr)
library(fields)
library(rootSolve)
library(patchwork)

source("R/0_helper-functions.R")
source("R/0_SIR-age-functions.R")

## some testing write SIR with time-varying beta (i.e., betahat)
ode_sir = function(t, y, parameters, bt) {
  with(as.list(c(y, parameters)), {
    # Define equations
    beta_seas = bt(t)
    dS = mu * (1-v) * N - beta_seas * S * I/N - S * mu
    dI = beta_seas * S * I/N - (mu + gamma) * I
    dR = mu * v * N + gamma * I - mu * R 
    res = c(dS, dI, dR)
    # Return list of gradients
    list(res)
  })
}


#### MODEL SETUP ---------------------------------------------------------------
source("R/2.1_setup-WAIFW.R") # this will add waifw to environment, which is list of WAIFW matrices to test

n_waifw = length(waifw)

compartments = c("S", "I", "R")

R0 <- 17
paras = c(mu = 1/80, N = 500000, beta1 = 0,
          gamma = 365/14, delta = 0, p = 0, phi = 1) # delta = 1e-4
paras["beta0"]  = (paras["gamma"] + paras["mu"])*R0

fert = rep(paras["mu"], length(age_classes))
mort = rep(paras["mu"], length(age_classes))

start_vax = 0.95

#### GET R0 VALUES FOR EACH WAIFW ----------------------------------------------
stable_age = findStableStruct(age_classes, mort, fert, 1/52)$stable.age

scalars = data.frame(waifw_id = 1:5, scalar = NA, diff = NA)
for(i in 1:length(waifw)){
  o = optimize(f = find_scalar, tol = 1e-8, interval = c(0, 100), R0 = R0,
               waifw = waifw[[i]], S = stable_age, beta0 = paras["beta0"], 
               gamma = paras["gamma"],  mu = paras["mu"], N = 1, age_classes = age_classes)
  print(get_Rt(waifw[[i]], stable_age, paras["beta0"]*o$minimum, paras["gamma"],  
               mu = paras["mu"], N = 1, age_classes = age_classes))
  scalars[i, 2:3] = c(o$minimum, o$objective)
}

# create a list of parameters for 
paras_all = lapply(1:5, function(i){paras_tmp = paras; paras_tmp["beta0"] = paras_tmp["beta0"]*scalars[i,2]; return(paras_tmp)})
paras_all_noimport = lapply(1:5, function(i){paras_tmp = paras; paras_tmp["beta0"] = paras_tmp["beta0"]*scalars[i,2]; return(paras_tmp)})

#### SET PRE-VACCINATION EQULIBRIUM --------------------------------------------
# pre-vax equilibrium = (1-vax_cov)*N(a) where N(a) is stable age distribution
vax_equilib_long = expand.grid(age = age_classes, 
                               start_vax = start_vax, 
                               variable = c("S", "I", "R")) %>%
  left_join(data.frame(age = age_classes, 
                       bin_width = bin_width, 
                       N = stable_age*paras["N"])) %>%
  mutate(value = ifelse(variable == "S", N*(1-start_vax), ifelse(variable == "I", 0, N*start_vax)))

#### SIMULATE FULL AGE STRUCTURED MODEL FOR EACH WAIFW -------------------------
#### INCLUDE BETA HAT ADJUSTMENT 
# with vaccination release
release_vax = 0.2
chosen_dt = 1/(52*7)
# show equilibrium values for different vax rates and waifw matrices
release_sim_df = vector("list", length(waifw))
# IC_release = vector("list", length(waifw))
for(i in 1:length(waifw)){
  print(paste0(i, "/", length(waifw)))
  IC_manual = vax_equilib_long %>% arrange(age)
  names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
  IC_manual = IC_manual$value
  names(IC_manual) = names_IC
  # try putting all I individuals in age 5
  I_indx = which(substr(names(IC_manual), 1, 1) == "I")
  IC_manual[which(names(IC_manual) == "I_5")] = 1
  IC_manual[which(names(IC_manual) == "R_5")] = IC_manual[which(names(IC_manual) == "R_5")] - 1
  IC_adj = vax_equilib_long %>% 
    summarize(value = sum(value), .by = c("variable")) %>% 
    pull(value) 
  names(IC_adj) = c("S_adj", "I_adj", "R_adj")
  IC_adj["I_adj"] = 1
  IC_adj["R_adj"] = IC_adj["R_adj"] - 1  
  release_sim_df[[i]] = run_ode(
    age_classes = age_classes, mort = mort,
    fert = fert, start_pop = paras["N"],
    compartments = compartments, dt = chosen_dt,
    vax_rates = c(release_vax), waifw = waifw[[i]],
    IC_manual = c(IC_manual, IC_adj), max_t = 10, params = paras_all[[i]], #, IC_manual = new_IC
    adjust_beta_flag = FALSE, 
    add_homogenous_adjusted_run = TRUE) %>%
    mutate(waifw_id = i)
}
beep()

release_sim_df_long = bind_rows(release_sim_df)

bh = release_sim_df_long %>%
  filter(variable %in% c("BH")) %>%
  select(-age) %>% 
  filter(!is.na(value)) %>%
  unique() %>%
  arrange(waifw_id, time)

# plot results
# beta hat
p1 = bh %>%
  mutate(waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))) %>%
  ggplot(aes(x = time, y = value/paras["gamma"], color = as.factor(waifw_id))) + 
  geom_hline(yintercept = paras["beta0"]/paras["gamma"], linetype = "dotted") +
  geom_line(linewidth = 0.5) + 
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since immunization decline") + 
  scale_y_continuous(name = expression(hat(beta)(t)/gamma)) +
  theme_bw(base_size = 7) +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())
# proportion infected
p2 = release_sim_df_long %>%
  filter(variable %in% c("I")) %>%
  summarize(value = sum(value), .by = c("time", "variable", "waifw_id")) %>%
  mutate(waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))) %>%
  ggplot(aes(x = time, y = value/paras["N"], color = as.factor(waifw_id))) + 
  geom_line(linewidth = 0.5) +
  geom_line(data = release_sim_df_long %>%
              filter(variable %in% c("I_adj"), !is.na(value)),
            linewidth = 0.5, color = "gray", linetype = "dashed") +
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since immunization decline") + 
  scale_y_continuous(name = "proportion infected") +
  theme_bw(base_size = 7) +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())
# difference between proportion infected in each case
p3 = release_sim_df_long %>%
  filter(variable %in% c("I")) %>%
  summarize(value = sum(value), .by = c("time", "variable", "waifw_id")) %>% 
  left_join(
    release_sim_df_long %>%
      filter(variable %in% c("I_adj"), !is.na(value)) %>%
      rename(value_adj = value) %>%
      select(-variable) 
  ) %>%
  mutate(I_diff = (value_adj-value)/value, 
         waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))) %>%
  ggplot(aes(x = time, y = I_diff, color = as.factor(waifw_id))) + 
  geom_line(linewidth = 0.5) +
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since immunization decline") + 
  scale_y_continuous(name = "percent difference in\nproportion infected") +
  theme_bw(base_size = 7) +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())
  
plot_grid(p1, p2, p3, ncol = 1, align = "v", axis = "lr", labels = c("A", "B", "C"), label_size = 8)
ggsave("figures/beta_hat_illustration.pdf", width = 6, height = 5)

#### SIMULATE STRUCTURED MODEL TO YIELD SIR DYNAMICS ---------------------------
# the equivalence ratio governs the transmission rate in a simple SIR model
# we expect to get the same dynamics as the flat WAIFW above

# use adjust_beta_flag = TRUE flag to implement time-varying beta_hat transmission rate
release_sim_df_adjbeta2 = vector("list", length(waifw))
for(i in 1:length(waifw)){
  print(paste0(i, "/", length(waifw)))
  IC_manual = vax_equilib_long %>% arrange(age)
  names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
  IC_manual = IC_manual$value
  names(IC_manual) = names_IC
  # try putting all I individuals in age 5
  I_indx = which(substr(names(IC_manual), 1, 1) == "I")
  IC_manual[which(names(IC_manual) == "I_5")] = 1
  IC_manual[which(names(IC_manual) == "R_5")] = IC_manual[which(names(IC_manual) == "R_5")] - 1
  release_sim_df_adjbeta2[[i]] = run_ode(
    age_classes = age_classes, mort = mort,
    fert = fert, start_pop = paras["N"],
    compartments = compartments, dt = chosen_dt,
    vax_rates = c(release_vax), waifw = waifw[[i]]*scalars[i, "scalar"],
    IC_manual = IC_manual, max_t = 10, params = paras_all[[1]], #, IC_manual = new_IC
    adjust_beta_flag = TRUE) %>%
    mutate(waifw_id = i)
}
beep()

release_sim_df_adjbeta2_long = bind_rows(release_sim_df_adjbeta2)

bh2 = release_sim_df_adjbeta2_long %>%
  filter(variable %in% c("BH")) %>%
  select(-age) %>% 
  filter(!is.na(value)) %>%
  unique() %>%
  arrange(waifw_id, time)

# plot results
p1 = bh2 %>%
  mutate(waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))) %>%
  ggplot(aes(x = time, y = paras_all[[1]]["beta0"]/value, color = as.factor(waifw_id))) + 
  geom_text(data = data.frame(y = c(0.95, 1/0.95), x = c(10, 10), vjust = c(1, 0),
                              waifw_id = factor(1, levels = c(1, 5, 4, 2, 3)),
                              lab = c("\nspeeding up\nwith age structure", "slowing down\nwith age structure\n")),
            aes(x = x, y = y, label = lab, vjust = vjust), 
            hjust = 1, color = "black", size = 1.9, alpha = 1, lineheight = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_line(linewidth = 0.5) + 
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since immunization decline") + 
  scale_y_log10(name = expression(k(t) ~ "(log scale)")) +
  theme_bw(base_size = 7) +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())
p2 = release_sim_df_long %>%
  filter(variable %in% c("I")) %>%
  summarize(value = sum(value), .by = c("time", "variable", "waifw_id")) %>%
  mutate(waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))) %>%
  ggplot(aes(x = time, y = value/paras["N"], color = as.factor(waifw_id))) + 
  geom_line(linewidth = 0.5, alpha = 0.4) +
  geom_line(data = release_sim_df_adjbeta2_long %>% filter(variable %in% c("I")) %>%
              summarize(value = sum(value), .by = c("time", "variable", "waifw_id")) %>%
              mutate(waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))),
            linewidth = 0.5, linetype = "dashed") +
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since immunization decline") + 
  scale_y_continuous(name = "proportion infected") +
  theme_bw(base_size = 7) +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())
p3 = release_sim_df_adjbeta2_long %>%
  filter(variable %in% c("I")) %>%
  summarize(value_adj = sum(value), .by = c("time", "variable", "waifw_id")) %>% 
  left_join(
    release_sim_df_adjbeta2_long %>%
      filter(variable == "I", waifw_id == 1) %>%
      summarize(value = sum(value), .by = c("time", "variable", "waifw_id")) %>% 
      select(-variable, -waifw_id) 
  ) %>%
  mutate(I_diff = (value_adj-value)/value, 
         waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))) %>%
  ggplot(aes(x = time, y = I_diff, color = as.factor(waifw_id))) + 
  geom_line(linewidth = 0.5) +
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since immunization decline") + 
  scale_y_continuous(name = "percent difference in\nproportion infected") +
  theme_bw(base_size = 7) +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())

cowplot::plot_grid(p1, p2, p3, ncol = 1, align = "v", axis = "lr", labels = c("A", "B", "C"), label_size = 8)
ggsave("figures/beta_hat_illustration_toSIR.pdf", width = 6, height = 5)

