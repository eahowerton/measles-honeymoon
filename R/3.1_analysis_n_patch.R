library(dplyr)
library(deSolve)
library(ggplot2)
library(cowplot)
library(reshape2)
library(beepr)

source("R/0_helper-functions.R")

# setup parameters for simulation
paras = c(mu = 1/80, N = 1, beta1 = 0, gamma = 365/14, delta = 0, R0 = 17)
paras["beta"] = paras["R0"] * (paras["gamma"] + paras["mu"])
pre_v = 0.95#c(0.94, 0.95, 0.96)

### BASIC SIR MODEL ------------------------------------------------------------

ode_sirs = function(t, y, parameters) {
  with(as.list(c(y, parameters)), {
    # Define equations
    dS = mu * (N) * (1 - v) - beta * S * I/N - mu * S
    dI = beta * S * I/N - (mu + gamma) * I
    dR = mu * (N) * v + gamma * I - mu * R
    res = c(dS, dI, dR)
    # Return list of gradients
    list(res)
  })
}

### RUN ANALYSIS ---------------------------------------------------------------

# set scenarios to test
start_vax = 0.95
new_vax = seq(0, 0.9, 0.01)
mu_to_test = paras["mu"]
c_tst = c(1, 0.99, 0.95, 0.9) 
step = 52
n_steps = 15*step

to_test = expand.grid(new_vax = new_vax,
                      mu = mu_to_test,
                      c = c_tst)

to_test_small = expand.grid(new_vax = new_vax,
                            mu = mu_to_test)

# get susceptible accumulation
susc_df = vector("list", nrow(to_test_small))
for(i in 1:nrow(to_test_small)){
  paras_tmp = paras
  paras_tmp["mu"] = to_test_small[i, "mu"]
  paras_tmp["beta"] = paras_tmp["R0"] * (paras_tmp["gamma"] + paras_tmp["mu"])
  paras_tmp["v"] = to_test_small[i, "new_vax"]
  inits_sirs = c(S = 1-start_vax, I = 0, R = start_vax)
  susc_df[[i]] = as.data.frame(deSolve::ode(inits_sirs, seq(0,n_steps/step,1/step), ode_sirs, paras_tmp))
}

susc_df = bind_rows(susc_df, .id = "s_id") %>%
  mutate(s_id = as.integer(s_id)) %>%
  left_join(to_test_small %>% mutate(s_id = seq_len(n()))) %>%
  select(-s_id, -I, -R)

# set number of patches (equivalent to 50%, 20% and 10% of patches drop)
tst_p = c(2, 5, 10)

susc_df2 = vector("list", length(tst_p))
for(i in 1:length(tst_p)){
  print(i)
  susc_df2[[i]] = to_test %>%
    full_join(susc_df, relationship = "many-to-many") %>%
    mutate(row_id = seq_len(n())) %>%
    mutate(Rt = get_Rt_npatch_dropvax(p = tst_p[i], S_drop = S, 
                                      start_vax = start_vax,
                                      beta0 =  paras["beta"],
                                      gamma = paras["gamma"],
                                      mu,
                                      c = c),
           .by = "row_id"
    ) %>%
    mutate(p = tst_p[i])
}
beepr::beep()

susc_df2 = bind_rows(susc_df2)

# now get honeymoon times
honeymoon_period = susc_df2 %>%
  filter(Rt > 1) %>%
  mutate(min_time = min(time), .by = c("new_vax", "mu", "c", "p")) %>%
  filter(time == min_time) %>%
  select(-min_time)

saveRDS(honeymoon_period, "data/output-data/honeymoon_period_npatch.rda")

# reconstruct model without patches
basic_sir_honeymoon = expand.grid(mu = paras["mu"],
                                  beta0 = paras["beta"],
                                  v = start_vax,
                                  gamma = paras["gamma"],
                                  new_v = seq(0, 0.90, 0.01)) %>%
  mutate(R0 = (beta0)/((gamma + mu)),
         herd_imm = 1-1/R0,
         S_star = ifelse(v > herd_imm, 1-v,
                         (mu + gamma)/beta0)) %>%
  mutate(time_to_cross = -1/mu * log(((1/R0) - (1-new_v))/(S_star - (1-new_v))))

ex_xs = c(0.1, 0.5, 0.75)
vert_seg = honeymoon_period %>%
  filter(new_vax %in% ex_xs, c %in% c(0.9, 1), p == 10) %>%
  mutate(c = paste0("c", c*100)) %>%
  select(c, p, time, new_vax) %>%
  tidyr::spread(c, time)
horiz_seg = honeymoon_period %>%
  filter(new_vax %in% ex_xs, c == 0.9, p == 10) %>%
  mutate(new_x = (exp(-mu*time)*(1-start_vax-1) - 1/17 + 1)/(1-exp(-mu*time)))
# (start_vax - exp(-mu*time)*(1/17-1))/(1-exp(-mu*time)))

p_labs =  paste0("coverage decrease in ", 1/tst_p*100, "% of patches")
names(p_labs) = tst_p

c_labs = paste("connectivity:", c_tst)
names(c_labs) = c_tst

p1 = ggplot(data = honeymoon_period %>% filter(p == 10), aes(x = new_vax, y = time, color = as.factor(1-c))) +
  geom_line(size = 0.5) +
  geom_segment(data = vert_seg, aes(x = new_vax, xend = new_vax, y = c100, yend = c90),
               color = "black", linetype = "dotted") +
  geom_segment(data = horiz_seg, aes(x = new_vax, xend = new_x, y = time, yend = time),
               color = "black", linetype = "dotted") +
  geom_point(data = honeymoon_period %>%
               filter(new_vax %in% ex_xs, c %in% c(0.9, 1), p == 10), show.legend = FALSE) +
  geom_point(data = horiz_seg, aes(x = new_x, y = time), color = "black", show.legend = FALSE) +
  # facet_wrap(vars(p), labeller = labeller(p = p_labs)) +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(4, "Oranges")[4:2]), name = "connectivity") +
  scale_x_continuous(breaks = seq(0, 9, 0.3), expand = c(0,0), limits = c(0,0.91),
                     label = scales::percent, name = "immunization after decline (in 10% of patches)") +
  scale_y_continuous(name = "theoretical honeymoon time\n(time to Re > 1)") +
  theme_bw(base_size = 7) +
  theme(legend.key.width = unit(0.4, "cm"),
        legend.position = "bottom", 
        strip.background = element_blank())
p1

p2 = ggplot(data = honeymoon_period %>% filter(c != 1), aes(x = new_vax, y = time, color = as.factor(p))) + 
  geom_line(size = 0.5) + 
  geom_line(data = honeymoon_period %>% filter(c == 1) %>% select(-c), color = "black", size = 0.5) + 
  facet_wrap(vars(c), ncol = 1, labeller = labeller(c = c_labs)) + 
  scale_color_manual(values = RColorBrewer::brewer.pal(4, "Purples")[4:2],
                     name = "patches\nin decline", labels = paste0(rev(100*1/tst_p), "%")) + 
  scale_x_continuous(breaks = seq(0, 0.9, 0.3), expand = c(0,0), limits = c(0,0.91),
                     label = scales::percent, name = "immunization after decline") +
  scale_y_continuous(name = "theoretical honeymoon time\n(time to Re > 1)") +
  theme_bw(base_size = 7) +
  theme(legend.key.width = unit(0.3, "cm"),
        legend.position = "bottom", 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())

plot_grid(p1, p2, rel_widths = c(0.6, 0.4), labels = c("A", "B"), label_size = 8)
ggsave("figures/connectivity_fig.pdf", width = 5, height = 3.5)

# values for text
honeymoon_period %>%
  filter(c != 1, p == 10) %>%
  left_join(basic_sir_honeymoon %>% rename(new_vax = new_v)) %>%
  mutate(diff = time - time_to_cross,
         pct_change = diff/time_to_cross) %>%
  summarize(mean = mean(pct_change),
            mn = min(pct_change),
            mx = max(pct_change), .by = c("c"))

