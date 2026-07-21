library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidytable)
library(cowplot)
library(dplyr)
library(beepr)
library(rootSolve)
library(patchwork)

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

fert = rep(paras["mu"], length(age_classes))
mort = rep(paras["mu"], length(age_classes))

start_vax = 0.95

#### GET R0 VALUES FOR EACH WAIFW ----------------------------------------------
stable_age = findStableStruct(age_classes, mort, fert, 1/52)$stable.age

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

#### SET PRE-VACCINATION EQULIBRIUM --------------------------------------------
# pre-vax equilibrium = (1-vax_cov)*N(a) where N(a) is stable age distribution
vax_equilib_long = expand.grid(age = age_classes, 
                               start_vax = start_vax, 
                               variable = c("S", "I", "R")) %>%
  left_join(data.frame(age = age_classes, 
                       bin_width = bin_width, 
                       N = stable_age*paras["N"])) %>%
  mutate(value = ifelse(variable == "S", N*(1-start_vax), ifelse(variable == "I", 0, N*start_vax)))

#### RELEASE VACCINATION WITH EACH WAIFW ---------------------------------------
tst = expand.grid(waifw_id = 1:5, 
                  release_vax = seq(0.2, 1, 0.2))
chosen_dt = 1/52
# show equilibrium values for different vax rates and waifw matrices
release_sim_df = vector("list", nrow(tst))
# IC_release = vector("list", length(waifw))
for(i in 1:nrow(tst)){
  print(paste0(i, "/", nrow(tst)))
  IC_manual = vax_equilib_long %>% arrange(age)
  names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
  IC_manual = IC_manual$value
  names(IC_manual) = names_IC
  # try putting all I individuals in age 5
  I_indx = which(substr(names(IC_manual), 1, 1) == "I")
  IC_manual[which(names(IC_manual) == "I_5")] = 1
  IC_manual[which(names(IC_manual) == "R_5")] = IC_manual[which(names(IC_manual) == "R_5")] - 1
  getIC = run_ode(
    age_classes = age_classes, mort = mort,
    fert = fert, start_pop = paras_jcm["N"],
    compartments = compartments, dt = chosen_dt,
    vax_rates = 0.01, waifw = waifw[[tst[i, "waifw_id"]]],
    IC_manual = IC_manual, max_t = 100, params = paras_all[[tst[i, "waifw_id"]]], #, IC_manual = new_IC
    adjust_beta_flag = FALSE) %>%
    mutate(waifw_id = tst[i, "waifw_id"], vax = tst[i, "release_vax"])
  IC_manual = getIC %>% filter(time == max(time), substr(variable,1,1) %in% c("S", "I", "R")) %>% arrange(age)
  names_IC = paste(IC_manual$variable, IC_manual$age, sep = "_")
  IC_manual = IC_manual$value
  names(IC_manual) = names_IC
  release_sim_df[[i]] = run_ode(
    age_classes = age_classes, mort = mort,
    fert = fert, start_pop = paras_jcm["N"],
    compartments = compartments, dt = chosen_dt,
    vax_rates = c(tst[i, "release_vax"]), waifw = waifw[[tst[i, "waifw_id"]]],
    IC_manual = IC_manual, max_t = 100, params = paras_all[[tst[i, "waifw_id"]]], #, IC_manual = new_IC
    adjust_beta_flag = FALSE) %>%
    mutate(waifw_id = tst[i, "waifw_id"], vax = tst[i, "release_vax"])
}

release_sim_df_long = bind_rows(release_sim_df)

ggplot(data = release_sim_df_long %>% 
         summarize(value = sum(value), .by = c("time", "variable", "waifw_id", "vax")) %>% 
         filter(variable == "I", time > 1, vax != 1), 
       aes(x = time, y = value, color = vax, id = )) + 
  geom_line() + 
  facet_grid(cols = vars(waifw_id))
  
release_sim_df_long %>% 
  filter(variable == "I", vax != 1) %>%
  mutate(year = floor(time)) %>%
  summarize(mean_age = sum(age*value)/sum(value), .by = c("year", "waifw_id", "vax")) %>% 
  mutate(waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))) %>%
  ggplot(aes(x = year, y = mean_age, color = as.factor(waifw_id),  group = as.factor(waifw_id))) + 
  geom_line() + 
  labs(x = "time since vax introduction", y = "mean age") +
  facet_wrap(vars(vax), labeller = label_both, nrow = 1) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.title = element_blank(), 
        legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave("/Users/emilyhowerton/Desktop/measles_mean_age_full.pdf", width = 8, height = 3)


release_sim_df_long %>% 
  filter(variable == "I", vax == 0.8, waifw_id %in% c(1, 2, 5)) %>%
  mutate(year = floor(time)) %>%
  summarize(mean_age = sum(age*value)/sum(value), .by = c("year", "waifw_id", "vax")) %>% 
  mutate(waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))) %>%
  ggplot(aes(x = year, y = mean_age, color = as.factor(waifw_id),  group = as.factor(waifw_id))) + 
  geom_line() + 
  labs(x = "years since vaccination introduction", y = "mean age of infection") +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  theme_bw() + 
  theme(legend.title = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank())
ggsave("/Users/emilyhowerton/Desktop/measles_mean_age.pdf", width = 5, height = 3)

