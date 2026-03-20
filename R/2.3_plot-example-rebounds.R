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
release_vax = 0.2
chosen_dt = 1/52
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
  release_sim_df[[i]] = run_ode(
    age_classes = age_classes, mort = mort,
    fert = fert, start_pop = paras_jcm["N"],
    compartments = compartments, dt = chosen_dt,
    vax_rates = c(release_vax), waifw = waifw[[i]],
    IC_manual = IC_manual, max_t = 10, params = paras_all[[i]], #, IC_manual = new_IC
    adjust_beta_flag = FALSE) %>%
    mutate(waifw_id = i)
}
beep()

release_sim_df_long = bind_rows(release_sim_df)

unity_beta_long = release_sim_df_long %>%
  filter(variable %in% c("BH")) %>%
  left_join(release_sim_df_long %>%
              filter(variable %in% c("BH"), waifw_id == 1) %>%
              select(time, variable, value) %>%
              rename(homog_value = value)) %>%
  left_join(scalars) %>%
  mutate(unity_beta = (scalar*value)/homog_value)

Rt_long = release_sim_df_long %>% 
  filter(variable == "S") %>%
  summarize(Rt = get_Rt(waifw = waifw[[waifw_id]], S = value, 
                        beta0 = paras_all[[waifw_id]]["beta0"], 
                        gamma = paras["gamma"], mu = paras["mu"], 
                        N = paras["N"]), .by = c("time", "waifw_id"))

max_contacts = lapply(waifw, melt) %>%
  bind_rows(.id = "waifw_id") %>%
  mutate(value_scaled = value/max(value), .by = c("waifw_id")) %>%
  mutate(waifw_id = as.integer(waifw_id)) %>%
  left_join(scalars %>% select(-diff)) %>% 
  mutate(value = value * scalar) %>% # rescale to have equivalent transmission rates
  mutate(m = max(value), .by = c("waifw_id")) %>%
  filter(value == m, waifw_id != 1) %>%
  mutate(waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))) %>%
  left_join(data.frame(Var1 = 1:length(age_classes), 
                       age1 = age_classes)) %>%
  left_join(data.frame(Var2 = 1:length(age_classes), 
                       age2 = age_classes))
  

#### PLOT ----------------------------------------------------------------------
pt0 = lapply(waifw_full, melt) %>%
  bind_rows(.id = "waifw_id") %>%
  mutate(value_scaled = value/max(value), .by = c("waifw_id")) %>%
  mutate(waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))) %>%
  left_join(data.frame(Var1 = 1:length(age_classes_full), 
                       age1 = age_classes_full)) %>%
  left_join(data.frame(Var2 = 1:length(age_classes_full), 
                       age2 = age_classes_full)) %>%
  ggplot(aes(x = age1, y = age2, fill = value_scaled)) + 
  geom_tile() + 
  geom_text(data = max_contacts, aes(label = round(m)), size = 1) +
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) +
  scale_fill_distiller(palette = "YlGnBu", name = "transmission\n(relative to max)") + 
  scale_x_continuous(expand = c(0,0),
                     name = "age (years)") +
  scale_y_continuous(expand = c(0,0),
                     name = "age (years)") +
  theme_gray(base_size = 7) +
  theme(#legend.position = "bottom", 
    legend.key.size = unit(0.34, "cm"),
        strip.background = element_blank())
ggsave("figures/waifws.pdf", pt0,  width = 8, height = 2.5)
pt1 = release_sim_df_long %>%
  filter(variable %in% c("I")) %>%
  summarize(value = sum(value), .by = c("time", "variable", "waifw_id")) %>%
  mutate(waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))) %>%
  ggplot(aes(x = time, y = value/paras["N"], color = as.factor(waifw_id))) + 
  geom_line(data = release_sim_df_long %>%
              filter(variable %in% c("I"), waifw_id == 1) %>%
              summarize(value = sum(value), .by = c("time", "variable", "waifw_id")) %>% select(-waifw_id), linewidth = 0.5, color = "black") +
  geom_line(linewidth = 0.5) + 
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since immunization decline") + 
  scale_y_continuous(name = "infections") +
  theme_bw(base_size = 7) +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())
pt2 = unity_beta_long %>%
  filter(variable == "BH") %>%
  mutate(waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))) %>%
  ggplot(aes(x = time, y = unity_beta, color = as.factor(waifw_id))) + 
  geom_text(data = data.frame(y = c(1e3, 1/1e3), x = c(10, 10), vjust = c(1, 0),
                              waifw_id = 1,
                              lab = c("\nspeeding up\nwith age structure", "slowing down\nwith age structure\n")),
            aes(x = x, y = y, label = lab, vjust = vjust), hjust = 1, color = "black", size = 1.9, alpha = 1) +
  geom_line(data = unity_beta_long %>%
              filter(variable == "BH", waifw_id == 1) %>% select(-waifw_id), linewidth = 0.5, color = "black") + 
  geom_line(aes(linetype = variable), linewidth = 0.5) + 
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  guides(color = FALSE) +
  coord_cartesian(ylim = c(1/1e3,1e3)) +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since immunization decline") + 
  scale_y_log10(breaks = c(1/rev(c(10, 1e2, 1e3)), c(1, 10, 1e2, 1e3)), 
                labels = c("1/1,000", "1/100", "1/10", "1", "10", "100", "1,000"), name = "equivalence ratio") +
  theme_bw(base_size = 7) +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())
pt3 = Rt_long %>%
  mutate(waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))) %>% 
  ggplot(aes(x = time, y = Rt, color = as.factor(waifw_id))) + 
  geom_hline(yintercept = 1, linetype = "dotted") + 
  geom_line(data = Rt_long %>% filter(waifw_id == 1) %>% select(-waifw_id), 
            linewidth = 0.5, color = "black") +
  geom_line(linewidth = 0.5) + 
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Set1")), labels = waifw_labs) +
  scale_x_continuous(breaks = seq(0,10,2), name = "years since immunization decline") + 
  scale_y_continuous(name = "Re") +
  theme_bw(base_size = 7) +
  theme(legend.title = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())
pt4 = release_sim_df_long %>% filter(variable %in% c("S", "I"), age <= 10) %>%
  mutate(bin_width = bin_width[as.factor(age)], 
         value = value/bin_width) %>%
  mutate(waifw_id = factor(waifw_id, levels = c(1, 5, 4, 2, 3))) %>% 
  dcast(time + age + waifw_id ~ variable, value.var = "value") %>%
  ggplot(aes(x = time, y = age)) + 
  geom_tile(aes(fill = S/paras["N"])) + 
  geom_point(data = release_sim_df_long %>% filter(time %in% seq(0, 10, 0.1), variable == "I", age <= 10) %>%
               mutate(bin_width = bin_width[as.factor(age)],
                      I = ifelse(value < 1, NA, value/bin_width)),
             aes(size = I), color = "white") +
  facet_grid(cols = vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) +
  guides(size = FALSE) + 
  scale_fill_viridis_c(option = "inferno", name = "% susceptible", labels = scales::percent) + 
  scale_size_continuous(range = c(0,1.25)) + 
  scale_x_continuous(expand = c(0,0), breaks = seq(0,10,2), name = "years since immunization decline") + 
  scale_y_continuous(expand = c(0,0), breaks = seq(0,10,2), name = "age (year)") + 
  theme_gray(base_size = 7) +
  theme(#legend.position = "none",
        legend.key.size = unit(0.35, "cm"),
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())

cowplot::plot_grid(pt0, pt3, pt1, pt4, pt2, ncol = 1, labels = LETTERS[1:5], align = "v", axis = "lr", label_size = 8)

ggsave("figures/release_examples_age.pdf", width = 6.5, height = 7)

# get values for text
# max unity beta in first period (not reliable for peak at age 10)
unity_beta_long %>% 
  left_join(data.frame(waifw_id = 1:5, 
                       approx_peak_time = c(3, 3, 3, 2,3))) %>%
  filter(time < approx_peak_time) %>% 
  summarize(m = 1/min(unity_beta, na.rm = TRUE), .by = waifw_id)


