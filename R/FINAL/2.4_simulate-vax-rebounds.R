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

start_vax = 0.95

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

#### SIMULATE REBOUND OUTBREAKS AT DIFFERENT STARTING TIMES --------------------
susc_after_release_long_full = readRDS("R/FINAL/DATA/susc_after_release_long.rds")

chosen_release_vax = 0.2
chosen_introduction_time = seq(1, 5, 0.5)
chosen_dt = 1/52
# show equilibrium values for different vax rates and waifw matrices
rebound_sim_df = vector("list", length(chosen_introduction_time))
# IC_release = vector("list", length(waifw))
for(j in 1:length(chosen_introduction_time)){
  tmp = vector("list", length(waifw))
  for(i in 1:length(waifw)){
    print(paste0(i, "/", length(waifw)))
    IC_manual = susc_after_release_long_full %>% filter(release_vax == chosen_release_vax, !(variable %in% c("BH", "BHs", "BHi", "BHb")), time == chosen_introduction_time[j])
    names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
    IC_manual = IC_manual$value
    names(IC_manual) = names_IC
    # try putting in 1 infected individual (say of age 5)
    I_indx = which(substr(names(IC_manual), 1, 1) == "I")
    tot_I = sum(IC_manual[I_indx])
    IC_manual[I_indx] = 0
    IC_manual[which(names(IC_manual) == "I_5")] = 1
    IC_manual[which(names(IC_manual) == "R_5")] = IC_manual[which(names(IC_manual) == "R_5")] - (tot_I - 1)
    tmp[[i]] = run_ode(
      age_classes = age_classes, mort = mort,
      fert = fert, start_pop = paras_jcm["N"],
      compartments = compartments, dt = chosen_dt,
      vax_change_times = c(0), vax_rates = c(chosen_release_vax), waifw = waifw[[i]],
      IC_type = "manual", IC_manual = IC_manual, max_t = 15, params = paras_all[[i]], #, IC_manual = new_IC
      adjust_beta_flag = FALSE, plot_flag = FALSE) %>%
      mutate(waifw_id = i, 
             introduction_time = chosen_introduction_time[j], 
             time = time + chosen_introduction_time[j] 
             )
  }
  rebound_sim_df[[j]] = bind_rows(tmp)
}
beep()

rebound_sim_df_long = bind_rows(rebound_sim_df)

#### NOW REINTRODUCE VACCINATION ONE YEAR LATER --------------------------------
chosen_release_vax = 0.2
chosen_introduction_time = seq(1, 5, 0.5)
revax_lag = 1
revax_level = 0.9
chosen_dt = 1/52
# show equilibrium values for different vax rates and waifw matrices
revax_sim_df = vector("list", length(chosen_introduction_time))
# IC_release = vector("list", length(waifw))
for(j in 1:length(chosen_introduction_time)){
  tmp = vector("list", length(waifw))
  for(i in 1:length(waifw)){
    print(paste0(i, "/", length(waifw)))
    IC_manual = rebound_sim_df_long %>% filter(waifw_id == i, !(variable %in% c("BH", "BHs", "BHi", "BHb")), 
                                               introduction_time == chosen_introduction_time[j], time == chosen_introduction_time[j] + revax_lag)
    names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
    IC_manual = IC_manual$value
    names(IC_manual) = names_IC
    tmp[[i]] = run_ode(
      age_classes = age_classes, mort = mort,
      fert = fert, start_pop = paras_jcm["N"],
      compartments = compartments, dt = chosen_dt,
      vax_change_times = c(0), vax_rates = c(revax_level), waifw = waifw[[i]],
      IC_type = "manual", IC_manual = IC_manual, max_t = 15, params = paras_all[[i]], #, IC_manual = new_IC
      adjust_beta_flag = FALSE, plot_flag = FALSE) %>%
      mutate(waifw_id = i, 
             introduction_time = chosen_introduction_time[j], 
             time = time + chosen_introduction_time[j] + revax_lag
      )
  }
  revax_sim_df[[j]] = bind_rows(tmp)
}
beep()


revax_sim_df_long = bind_rows(revax_sim_df)

revax_sim_plot = revax_sim_df_long %>%
  filter(variable %in% c("I")) %>%
  summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "introduction_time"))
  
rebound_sim_df_long %>%
  filter(variable %in% c("I"), introduction_time == floor(introduction_time)) %>%
  summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "introduction_time")) %>%
  ggplot(aes(x = time, y = value, color = as.factor(waifw_id))) + 
  geom_vline(aes(xintercept = introduction_time)) +
  geom_vline(aes(xintercept = introduction_time + revax_lag), linetype = "dotted") +
  geom_line(data = revax_sim_plot %>% filter(introduction_time == floor(introduction_time)), alpha = 0.5, linewidth = 0.8) +
  geom_line(linewidth = 0.8) + 
  facet_grid(cols = vars(waifw_id), rows = vars(introduction_time), labeller = labeller(waifw_id = waifw_labs), scales = "free") + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since release", limits = c(0, 10)) + 
  scale_y_continuous(name = "infections") +
  theme_bw() + 
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())
# ggsave("R/FINAL/figures/introduction_timing.pdf", width = 8, height = 6)


rebound_sim_df_long %>%
  filter(variable %in% c("I"), introduction_time == 3) %>%
  summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "introduction_time")) %>%
  ggplot(aes(x = time, y = value, color = as.factor(waifw_id))) + 
  geom_vline(aes(xintercept = introduction_time)) +
  geom_vline(aes(xintercept = introduction_time + revax_lag), linetype = "dotted") +
  # geom_line(data = rebound_sim_df_long %>%
  #             filter(variable %in% c("I"), waifw_id == 1, introduction_time == 3) %>%
  #             summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "introduction_time")) %>% select(-waifw_id), linewidth = 0.8, color = "black") +
  geom_line(data = revax_sim_plot %>% filter(introduction_time == 3), alpha = 0.5, linewidth = 0.8) +
  geom_line(linewidth = 0.8) + 
  facet_grid(cols = vars(waifw_id), rows = vars(introduction_time), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since release", limits = c(0, 10)) + 
  scale_y_continuous(name = "infections") +
  theme_bw() + 
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())

# total infections in 15 year period after introduction
p1 = revax_sim_plot %>% 
  filter(variable == "I") %>%
  summarize(total_I = sum(value), 
            nyears = max(time) - min(time), .by = c("waifw_id", "introduction_time")) %>%
  # need to add first year of cases from non-vax period
  left_join(
    rebound_sim_df_long %>%
      filter(variable %in% c("I"), time <= introduction_time + 1) %>%
      summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "introduction_time")) %>%
      filter(variable == "I") %>%
      summarize(total_I_prevax_period = sum(value), .by = c("waifw_id", "introduction_time"))
  ) %>% 
  mutate(total_I = total_I + total_I_prevax_period) %>%
  select(-total_I_prevax_period) %>%
  mutate(revax = TRUE) %>%
  bind_rows(rebound_sim_df_long %>%
              filter(variable %in% c("I")) %>%
              summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "introduction_time")) %>%
              filter(variable == "I") %>%
              summarize(total_I = sum(value),
                        nyears = max(time) - min(time), .by = c("waifw_id", "introduction_time")) %>%
              mutate(revax = FALSE)) %>%
  ggplot(aes(x = introduction_time, y = total_I, color = as.factor(waifw_id), linetype = revax)) + 
  geom_point() + 
  geom_line(alpha = 0.7) +
  facet_grid(rows = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  guides(color = "none") + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())

p2 = revax_sim_plot %>% 
  filter(variable == "I") %>%
  summarize(total_I = cumsum(value), .by = c("waifw_id", "introduction_time", "time")) %>%
  mutate(revax = TRUE) %>%
  bind_rows(rebound_sim_df_long %>%
              filter(variable %in% c("I")) %>%
              summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "introduction_time")) %>%
              filter(variable == "I") %>%
              summarize(total_I = cumsum(value), .by = c("waifw_id", "introduction_time", "time")) %>%
              mutate(revax = FALSE)) %>%
  filter(introduction_time %in% c(1,3,5)) %>%
  ggplot(aes(x = time, y = total_I, color = as.factor(waifw_id), alpha = revax)) + 
  geom_vline(data = data.frame(introduction_time = c(1,3,5)), aes(xintercept = introduction_time)) +
  geom_vline(data = data.frame(introduction_time = c(1,3,5), 
                               vax_start_time = c(1,3,5)+1), aes(xintercept = vax_start_time), linetype = "dashed") +  # vax introduction
  geom_line() +
  facet_grid(rows = vars(waifw_id), cols = vars(introduction_time), labeller = labeller(waifw_id = waifw_labs), scales = "free") + 
  guides(color = "none") + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_alpha_manual(values = c(0.5, 1)) +
  theme_bw() + 
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())

plot_grid(p2, p1, rel_widths = c(0.7, 0.3), align = "h", axis = "tb")
ggsave("R/FINAL/figures/averted_cases.pdf", width = 12, height = 7)


# total infections in 5 years after introduction + 1 year revax lag
revax_sim_plot %>% 
  filter(variable == "I", time < introduction_time + revax_lag + 5) %>%
  summarize(total_I = sum(value), 
            nyears = max(time) - min(time), .by = c("waifw_id", "introduction_time")) %>%
  mutate(revax = TRUE) %>%
  bind_rows(rebound_sim_df_long %>%
              filter(variable %in% c("I"), time < introduction_time + revax_lag + 5) %>%
              summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "introduction_time")) %>% 
              filter(variable == "I") %>%
              summarize(total_I = sum(value), 
                        nyears = max(time) - min(time), .by = c("waifw_id", "introduction_time")) %>%
              mutate(revax = FALSE)) %>%
  ggplot(aes(x = introduction_time, y = total_I, color = as.factor(waifw_id), linetype = revax)) + 
  geom_point() + 
  geom_line(alpha = 0.7) +
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  guides(color = "none") + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())

# averted 
revax_sim_plot %>% 
  filter(variable == "I", time < introduction_time + revax_lag + 5) %>%
  summarize(total_I = sum(value), 
            nyears = max(time) - min(time), .by = c("waifw_id", "introduction_time")) %>%
  mutate(revax = TRUE) %>%
  bind_rows(rebound_sim_df_long %>%
              filter(variable %in% c("I"), time < introduction_time + revax_lag + 5) %>%
              summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "introduction_time")) %>% 
              filter(variable == "I") %>%
              summarize(total_I = sum(value), 
                        nyears = max(time) - min(time), .by = c("waifw_id", "introduction_time")) %>%
              mutate(revax = FALSE)) %>%
  mutate(revax = paste0("revax_", revax)) %>%
  dcast(waifw_id + introduction_time ~ revax, value.var = "total_I") %>%
  ggplot(aes(x = introduction_time, y = (revax_TRUE - revax_FALSE)/5, color = as.factor(waifw_id))) + 
  geom_point() + 
  geom_line(alpha = 0.7) +
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  guides(color = "none") + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())

rebound_sim_df_long %>%
  filter(variable %in% c("I"), introduction_time == floor(introduction_time), time > introduction_time + revax_lag) %>%
  summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "introduction_time")) %>%
  mutate(cum_value = cumsum(value), .by = c("variable", "waifw_id", "introduction_time")) %>%
  ggplot(aes(x = time, y = cum_value, color = as.factor(waifw_id))) + 
  geom_vline(aes(xintercept = introduction_time)) +
  geom_vline(aes(xintercept = introduction_time + revax_lag), linetype = "dotted") +
  geom_line(data = revax_sim_df_long %>%
              filter(variable %in% c("I")) %>%
              summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "introduction_time")) %>%
              mutate(cum_value = cumsum(value), .by = c("variable", "waifw_id", "introduction_time")) %>%
              filter(introduction_time == floor(introduction_time)), alpha = 0.5, linewidth = 0.8) +
  geom_line(linewidth = 0.8) + 
  facet_grid(cols = vars(waifw_id), rows = vars(introduction_time), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since release", limits = c(0, 10)) + 
  scale_y_continuous(name = "infections") +
  theme_bw() + 
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())
  
rebound_sim_df_long %>%
  filter(variable %in% c("I"), introduction_time %in% c(1,3,5), time > introduction_time + revax_lag) %>%
  summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "introduction_time")) %>%
  ggplot(aes(x = time, y = value, color = as.factor(waifw_id))) + 
  geom_vline(aes(xintercept = introduction_time)) +
  geom_vline(aes(xintercept = introduction_time + revax_lag), linetype = "dotted") +
  geom_line(data = revax_sim_df_long %>%
              filter(variable %in% c("I")) %>%
              summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "introduction_time")) %>%
              filter(introduction_time %in% c(1,3,5)), linewidth = 0.8) +
  geom_line(data = rebound_sim_df_long %>%
              filter(variable %in% c("I"), introduction_time %in% c(1,3,5), time >= introduction_time, time <= introduction_time + revax_lag) %>%
              summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "introduction_time"))  %>%
              filter(introduction_time %in% c(1,3,5)), linewidth = 0.8) +
  geom_line(linewidth = 0.8, alpha = 0.5) + 
  facet_grid(cols = vars(waifw_id), rows = vars(introduction_time), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since release", limits = c(0, 10)) + 
  scale_y_continuous(name = "infections") +
  theme_bw() + 
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())
ggsave("R/FINAL/figures/introduction_timing.pdf", width = 8, height = 4)

