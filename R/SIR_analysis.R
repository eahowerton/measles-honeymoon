## try to recreate Bryan's SIR results
## start with Bjornstad example (section 6.3)

library(dplyr)
library(deSolve)
library(ggplot2)
library(reshape2)

seirmod2 = function(t, y, parameters, vax_pct) {
  S = y[1]
  E = y[2]
  I = y[3]
  R = y[4]
  parameters = c(parameters, v = vax_pct(t))
  with(as.list(parameters), {
    beta = beta0 * (1 + beta1 * cos(2 * pi * t))
    dS = mu * (N - S) * (1 - v) - beta * S * I/N - delta * S # additional importations
    dE = beta * S * I/N - (mu + sigma) * E + delta * S
    dI = sigma * E - (mu + gamma) * I
    dR = mu * (N - S) * v + gamma * I - mu * R
    res = c(dS, dE, dI, dR)
    list(res)
  })
}

# times = c(seq(0, 80-1/365, by = 1/365), seq(82, 85, by = 1/(365*5)), seq(85+1/365, 100, by=1/365))
times = seq(0, 100, by = 1/365)
no_vax = approxfun(times, rep(0, length(times)))
paras = c(mu = 1/50, N = 1, beta0 = 1000, beta1 = 0.2,
          sigma = 365/8, gamma = 365/5, vax_pct = 0, delta = 0)
start = c(S = 0.06, E = 0, I = 0.001, R = 0.939)
results = as.data.frame(ode(start, times, seirmod2, paras, vax_pct = no_vax))

ggplot(data = results %>% filter(time > 80), aes(x = time, y = I)) + 
  geom_line()

# quick bifurcation diagram
beta1_vec = seq(0, 0.25, length = 101)

bifur_reslts = list()
# Loop over beta1’s
for (i in 1:length(beta1_vec)) {
  paras_tmp = paras
  paras_tmp["beta1"] = beta1_vec[i]
  bifur_reslts[[i]] = as.data.frame(ode(start, times, seirmod2, paras_tmp, vax_pct = no_vax)) %>%
    mutate(beta1 = beta1_vec[i])
}

# based on bifurcation results, choose three examples
annual_beta1 = 0.1
bienn1_beta1 = 0.16
bienn2_beta1 = 0.2

examp_beta1 = c(annual_beta1, bienn1_beta1, bienn2_beta1)

p1 = bind_rows(bifur_reslts) %>%
  filter(time%%1 == 0.25, time > 80) %>%
  ggplot(aes(x = beta1, y = I)) + 
  geom_point(alpha = 0.2) + 
  geom_vline(data = data.frame(beta1 = examp_beta1), 
             aes(xintercept = beta1), linetype = "dotted", color = 'red') + 
  theme_bw()
p2 = bind_rows(bifur_reslts) %>%
  filter(round(beta1, 5) %in% examp_beta1, time > 90) %>%
  ggplot(aes(x = time, y = I)) + 
  geom_line() + 
  facet_wrap(vars(beta1), ncol = 1) + 
  theme_bw()
p3 = bind_rows(bifur_reslts) %>%
  filter(round(beta1, 5) %in% examp_beta1, time > 94) %>%
  ggplot(aes(x = S, y = I)) + 
  geom_path() + 
  facet_wrap(vars(beta1), ncol = 1) + 
  theme_bw()
cowplot::plot_grid(p1, p2, p3, nrow = 1)

## now implement vaccination during this period
bind_rows(bifur_reslts) %>%
  filter(beta1 %in% examp_beta1, time > 81, time < 82) %>%
  mutate(max_I = max(I), .by = c("beta1")) %>%
  filter(I == max_I)

vax_on = 0.9
vax_release = 0.4
nyear_vax = 3
time_on_peak = 80.25
time_on_trough = 81.25

# vax_peakstart = approxfun(times, c(rep(0, time_on_peak*120+1), 
#                                    rep(vax_on, nyear_vax*120), 
#                                    rep(vax_release, length(times) - time_on_peak*120 - 1 - nyear_vax*120)))
# vax_troughstart = approxfun(times, c(rep(0, time_on_trough*120+1), 
#                                      rep(vax_on, nyear_vax*120), 
#                                      rep(vax_release, length(times) - time_on_trough*120 - 1 - nyear_vax*120)))

vax_peakstart = approxfun(c(0, time_on_peak-1/365, time_on_peak, time_on_peak + nyear_vax, time_on_peak + nyear_vax + 1/365,  max(times)), 
                          c(0, 0, vax_on, vax_on, vax_release, vax_release))
vax_troughstart = approxfun(c(0, time_on_trough-1/365, time_on_trough, time_on_trough + nyear_vax, time_on_trough + nyear_vax + 1/365, max(times)), 
                            c(0, 0, vax_on, vax_on, vax_release, vax_release))


vax_peak_results = list()
vax_trough_results = list()
for(i in 1:length(examp_beta1)){
  paras_tmp = paras
  paras_tmp["beta1"] = examp_beta1[i]
  vax_peak_results[[i]] = as.data.frame(ode(start, times, seirmod2, paras_tmp, vax_pct = vax_peakstart, rtol = 1e-10, atol = 1e-10)) %>%
    mutate(beta1 = examp_beta1[i], 
           # Re = c(0, diff(Re)), 
           Re = paras_tmp["beta0"] * (1 + paras_tmp["beta1"] * cos(2 * pi * time)) * paras_tmp["sigma"] * S / ((paras_tmp["gamma"] + paras_tmp["mu"]) * (paras_tmp["sigma"] + paras_tmp["mu"])  * paras_tmp["N"]))
  vax_trough_results[[i]] = as.data.frame(ode(start, times, seirmod2, paras_tmp, vax_pct = vax_troughstart)) %>%
    mutate(beta1 = examp_beta1[i], 
           Re = paras_tmp["beta0"] * (1 + paras_tmp["beta1"] * cos(2 * pi * time)) * paras_tmp["sigma"] * S / ((paras_tmp["gamma"] + paras_tmp["mu"]) * (paras_tmp["sigma"] + paras_tmp["mu"])  * paras_tmp["N"]))
}

### now what if there's a slow decline in vax rates rather than a sudden jump
nyear_return = 5
vax_peakstart_grad = approxfun(c(0, time_on_peak-1/120, time_on_peak, time_on_peak + nyear_vax, time_on_peak + nyear_vax + nyear_return, max(times)), 
                          c(0, 0, vax_on, vax_on, vax_release, vax_release))
vax_troughstart_grad = approxfun(c(0, time_on_trough-1/120, time_on_trough, time_on_trough + nyear_vax, time_on_trough + nyear_vax + nyear_return, max(times)), 
                          c(0, 0, vax_on, vax_on, vax_release, vax_release))

vax_peak_results_grad = list()
vax_trough_results_grad = list()
for(i in 1:length(examp_beta1)){
  paras_tmp = paras
  paras_tmp["beta1"] = examp_beta1[i]
  vax_peak_results_grad[[i]] = as.data.frame(ode(start, times, seirmod2, paras_tmp, vax_pct = vax_peakstart_grad)) %>%
    mutate(beta1 = examp_beta1[i], 
           Re = paras_tmp["beta0"] * (1 + paras_tmp["beta1"] * cos(2 * pi * time)) * paras_tmp["sigma"] * S / ((paras_tmp["gamma"] + paras_tmp["mu"]) * (paras_tmp["sigma"] + paras_tmp["mu"])  * paras_tmp["N"]))
  vax_trough_results_grad[[i]] = as.data.frame(ode(start, times, seirmod2, paras_tmp, vax_pct = vax_troughstart_grad)) %>%
    mutate(beta1 = examp_beta1[i], 
           Re = paras_tmp["beta0"] * (1 + paras_tmp["beta1"] * cos(2 * pi * time)) * paras_tmp["sigma"] * S / ((paras_tmp["gamma"] + paras_tmp["mu"]) * (paras_tmp["sigma"] + paras_tmp["mu"])  * paras_tmp["N"]))
}


rslts_all = bind_rows(vax_peak_results_grad) %>%
  mutate(vax_start = "peak", honeymoon_end = "gradual") %>%
  bind_rows(bind_rows(vax_trough_results_grad) %>% mutate(vax_start = "trough", honeymoon_end = "gradual")) %>%
  bind_rows(bind_rows(vax_trough_results) %>% mutate(vax_start = "trough", honeymoon_end = "jump")) %>%
  bind_rows(bind_rows(vax_peak_results) %>% mutate(vax_start = "peak", honeymoon_end = "jump")) %>%
  mutate(stage = ifelse(time < ifelse(vax_start == "peak", time_on_peak, time_on_trough), "pre vax", 
                        ifelse(time < ifelse(vax_start == "peak", time_on_peak, time_on_trough) + nyear_vax, "high vax", "decline vax")), 
         time = ifelse(vax_start == "trough", time - 1, time)
         ) 


p1 = data.frame(
  time = c(times, times), 
  peak = c(vax_peakstart_grad(times), vax_peakstart(times)),
  trough = c(vax_troughstart_grad(times), vax_troughstart(times)), 
  honeymoon_end = c(rep("graudual", length(times)), rep("jump", length(times)))
) %>%
  reshape2::melt(c("time", "honeymoon_end")) %>%
  filter(time > 77) %>%
  mutate(time = ifelse(variable == "trough", time - 1, time)) %>%
  ggplot(aes(x = time, y = value)) + 
  geom_line() + 
  facet_grid(cols = vars(honeymoon_end)) + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  scale_y_continuous(limits = c(0,1), 
                     name = "vax coverage") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = "none")

p2 = rslts_all %>%
  filter(time > 77) %>%
  ggplot() + 
  geom_rect(data = data.frame(t = c(time_on_peak),
                              tend = c(time_on_peak) + nyear_vax,
                              vax_start = c("peak")), 
            aes(xmin = t, xmax = tend, ymin = 0, ymax = Inf), alpha = 0.2) +
  geom_vline(data = data.frame(t = seq(80, 100, 5)), aes(xintercept = t), linetype = "dotted") + 
  geom_line(aes(x = time, y = I, color = stage, linetype = vax_start)) + 
  facet_grid(rows = vars(beta1), cols = vars(honeymoon_end)) + 
  scale_color_manual(values = c("blue", "red", "black")) + 
  scale_linetype_manual(values = c("longdash", "solid")) + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank())

cowplot::plot_grid(p1, p2, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.3, 0.7))

rslts_all %>%
  filter(time > 77, time < 85) %>%
  ggplot(aes(x = log(S), y = log(I), color = stage, linetype = vax_start)) + 
  geom_path() + 
  ggtitle("without importation") + 
  facet_grid(rows = vars(beta1), cols = vars(honeymoon_end)) + 
  scale_color_manual(values = c("blue", "red", "black")) + 
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank())


rslts_all %>%
  filter(time > 77) %>%
  ggplot(aes(x = time, y = log(Re), color = stage, linetype = vax_start)) + 
  geom_hline(yintercept = 0) + 
  # geom_rect() +
  geom_path() + 
  facet_grid(rows = vars(beta1), cols = vars(honeymoon_end)) + 
  ggtitle("without importation") + 
  scale_color_manual(values = c("blue", "red", "black")) + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank())

rslts_all %>%
  filter(time > 77, time < 85) %>%
  ggplot(aes(x = log(Re), y = log(I), color = stage, linetype = vax_start)) + 
  # geom_rect() +
  geom_path() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  facet_grid(rows = vars(beta1), cols = vars(honeymoon_end)) + 
  scale_color_manual(values = c("blue", "red", "black")) + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank())

rslts_all %>%
  filter(time > 77, time < 85, honeymoon_end == "jump") %>%
  mutate(stage = ifelse(time < ifelse(vax_start == "peak", time_on_peak, time_on_trough), "pre vax", 
                        ifelse(time < ifelse(vax_start == "peak", time_on_peak, time_on_trough) + nyear_vax, "high vax", "decline vax"))) %>%
  ggplot(aes(x = log(Re), y = log(I), color = stage, linetype = vax_start)) + 
  # geom_rect() +
  geom_path() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  facet_grid(rows = vars(beta1)) + 
  scale_color_manual(values = c("blue", "red", "black")) + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank())


# ## try some 3D plots
library(plotly)

p <- plot_ly(data = bind_rows(vax_peak_results) %>% 
               mutate(vax_start = "peak", honeymoon_end = "jump") %>%
               bind_rows(bind_rows(vax_trough_results) %>% mutate(vax_start = "trough", honeymoon_end = "jump")) %>%
               mutate(stage = ifelse(time < ifelse(vax_start == "peak", time_on_peak, time_on_trough), "pre vax", 
                                     ifelse(time < ifelse(vax_start == "peak", time_on_peak, time_on_trough) + nyear_vax, 
                                            "high vax", "decline vax")), 
                      time = ifelse(vax_start == "trough", time - 1, time)) %>%
               filter(time > 75, beta1 == 0.2),
             x = ~time, y = ~S, z = ~I, color = ~stage, linetype = ~vax_start, type = 'scatter3d', mode = 'lines', 
             linetypes = c(gradual = "dotted", jump = "solid"), 
             colors = c(`pre vax` = "black", `high vax` = "red", `decline vax` = "blue")) 
p 

p2 <- plot_ly(data = bind_rows(vax_trough_results_grad) %>% 
               mutate(vax_start = "trough", honeymoon_end = "gradual") %>%
               bind_rows(bind_rows(vax_trough_results) %>% mutate(vax_start = "trough", honeymoon_end = "jump")) %>%
               mutate(stage = ifelse(time < ifelse(vax_start == "peak", time_on_peak, time_on_trough), "pre vax", 
                                     ifelse(time < ifelse(vax_start == "peak", time_on_peak, time_on_trough) + nyear_vax, 
                                            "high vax", "decline vax"))) %>%
               filter(time > 75, time < 90, beta1 == 0.2),
             x = ~time, y = ~log(Re), z = ~I, color = ~stage, linetype = ~honeymoon_end, type = 'scatter3d', mode = 'lines', 
             linetypes = c(gradual = "dotted", jump = "solid"), 
             colors = c(`pre vax` = "black", `high vax` = "red", `decline vax` = "blue")) 
p2

## okay, so the hypothesis is that the driver of transient differneces is 
## I at the point where Re = 1, so let's try to set I = 1 (or some fixed value)
## in all cases
## another idea: add some importation to the model and see if anything changes

new_ICs = rslts_all %>% 
  filter(time > time_on_peak + nyear_vax, 
         Re > 1) %>%
  mutate(min_time = min(time), .by = c("honeymoon_end", "beta1", "vax_start")) %>%
  filter(time == min_time) %>%
  mutate(new_I = 1e-6, 
         R = I - new_I + R, 
         I = new_I,
         N = S + E + I + R)

max_t = new_ICs %>% pull(time) %>% max()
  
# now rerurn
vax_peak_results_updateIC = list()
vax_trough_results_updateIC = list()
vax_peak_results_grad_updateIC = list()
vax_trough_results_grad_updateIC = list()
for(i in 1:length(examp_beta1)){
  paras_tmp = paras
  paras_tmp["beta1"] = examp_beta1[i]
  # peak start jump decline
  vax_peak_results_updateIC[[i]] = as.data.frame(
    ode(new_ICs %>% filter(vax_start == "peak", honeymoon_end == "jump", beta1 == examp_beta1[i]) %>%
          select(S, E, I, R) %>% unlist(), 
        times[which(times > max_t)], seirmod2, paras_tmp, vax_pct = no_vax)) %>%
    mutate(beta1 = examp_beta1[i], 
           Re = paras_tmp["beta0"] * (1 + paras_tmp["beta1"] * cos(2 * pi * time)) * 
             paras_tmp["sigma"] * S / ((paras_tmp["gamma"] + paras_tmp["mu"]) * (paras_tmp["sigma"] + paras_tmp["mu"])
                                       * paras_tmp["N"]))
  # trough start jump decline
  vax_trough_results_updateIC[[i]] = as.data.frame(
    ode(new_ICs %>% filter(vax_start == "trough", honeymoon_end == "jump", beta1 == examp_beta1[i]) %>%
          select(S, E, I, R) %>% unlist(), 
        times[which(times > max_t)], seirmod2, paras_tmp, vax_pct = no_vax)) %>%
    mutate(beta1 = examp_beta1[i], 
           Re = paras_tmp["beta0"] * (1 + paras_tmp["beta1"] * cos(2 * pi * time)) * 
             paras_tmp["sigma"] * S / ((paras_tmp["gamma"] + paras_tmp["mu"]) * (paras_tmp["sigma"] + paras_tmp["mu"])  * 
                                         paras_tmp["N"]))
  # peak start gradual decline
  vax_peak_results_grad_updateIC[[i]] = as.data.frame(
    ode(new_ICs %>% filter(vax_start == "peak", honeymoon_end == "gradual", beta1 == examp_beta1[i]) %>%
          select(S, E, I, R) %>% unlist(), 
        times[which(times > max_t)], seirmod2, paras_tmp, vax_pct = no_vax)) %>%
    mutate(beta1 = examp_beta1[i], 
           Re = paras_tmp["beta0"] * (1 + paras_tmp["beta1"] * cos(2 * pi * time)) * 
             paras_tmp["sigma"] * S / ((paras_tmp["gamma"] + paras_tmp["mu"]) * (paras_tmp["sigma"] + paras_tmp["mu"])
                                       * paras_tmp["N"]))
  # trough start gradual decline
  vax_trough_results_grad_updateIC[[i]] = as.data.frame(
    ode(new_ICs %>% filter(vax_start == "trough", honeymoon_end == "gradual", beta1 == examp_beta1[i]) %>%
          select(S, E, I, R) %>% unlist(), 
        times[which(times > max_t)], seirmod2, paras_tmp, vax_pct = no_vax)) %>%
    mutate(beta1 = examp_beta1[i], 
           Re = paras_tmp["beta0"] * (1 + paras_tmp["beta1"] * cos(2 * pi * time)) * 
             paras_tmp["sigma"] * S / ((paras_tmp["gamma"] + paras_tmp["mu"]) * (paras_tmp["sigma"] + paras_tmp["mu"])  * 
                                         paras_tmp["N"]))
}

rslts_all_update_IC <- bind_rows(vax_peak_results_grad_updateIC) %>%
  mutate(vax_start = "peak", honeymoon_end = "gradual") %>%
  bind_rows(bind_rows(vax_trough_results_grad_updateIC) %>% mutate(vax_start = "trough", honeymoon_end = "gradual")) %>%
  bind_rows(bind_rows(vax_trough_results_updateIC) %>% mutate(vax_start = "trough", honeymoon_end = "jump")) %>%
  bind_rows(bind_rows(vax_peak_results_updateIC) %>% mutate(vax_start = "peak", honeymoon_end = "jump")) %>%
  mutate(stage = ifelse(time < ifelse(vax_start == "peak", time_on_peak, time_on_trough), "pre vax", 
                        ifelse(time < ifelse(vax_start == "peak", time_on_peak, time_on_trough) + nyear_vax, "high vax", "decline vax")))

rslts_all_update_IC %>%
  filter(time > 77) %>%
  ggplot() + 
  geom_rect(data = data.frame(t = c(time_on_peak, time_on_trough),
                              tend = c(time_on_peak, time_on_trough) + nyear_vax,
                              vax_start = c("peak", "trough")), 
            aes(xmin = t, xmax = tend, ymin = 0, ymax = Inf), alpha = 0.2) +
  geom_line(aes(x = time, y = I, color = stage, linetype = vax_start)) + 
  geom_text(data = new_ICs %>% dcast(beta1 + honeymoon_end ~ vax_start, value.var = "Re") %>% mutate(text = paste0("peak Re: ", round(peak, 5), "\ntrough Re: ", round(trough, 5))), 
            aes(x = Inf, y = Inf, label = text), hjust = 1, vjust = 1)+ 
  facet_grid(rows = vars(beta1), cols = vars(honeymoon_end)) + 
  scale_color_manual(values = c("blue", "red", "black")) + 
  scale_linetype_manual(values = c("longdash", "solid")) + 
  scale_y_continuous(limits = c(0, 0.005)) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank())
  
#### ADD SMALL AMOUNT OF IMPORTATION -------------------------------------------
paras_importation = c(mu = 1/50, N = 1, beta0 = 1000, beta1 = 0.2,
          sigma = 365/8, gamma = 365/5, vax_pct = 0, delta = 1e-4) # delta chosen arbitrarily for now

results_importation = as.data.frame(ode(start, times, seirmod2, paras_importation, vax_pct = no_vax))

ggplot(data = results %>% filter(time > 95), aes(x = time, y = I)) + 
  geom_line() + 
  geom_line(data = results_importation %>% filter(time > 95), color = "red")

# quick bifurcation diagram again
bifur_reslts_import = list()
# Loop over beta1’s
for (i in 1:length(beta1_vec)) {
  paras_tmp = paras_importation
  paras_tmp["beta1"] = beta1_vec[i]
  bifur_reslts_import[[i]] = as.data.frame(ode(start, times, seirmod2, paras_tmp, vax_pct = no_vax)) %>%
    mutate(beta1 = beta1_vec[i])
}

p1 = bind_rows(bifur_reslts_import) %>%
  filter(time%%1 == 0, time > 80, !is.na(I)) %>%
  ggplot(aes(x = beta1, y = I)) + 
  geom_point(alpha = 0.2) + 
  geom_vline(data = data.frame(beta1 = examp_beta1), 
             aes(xintercept = beta1), linetype = "dotted", color = 'red') + 
  theme_bw()
p2 = bind_rows(bifur_reslts_import) %>%
  filter(round(beta1, 5) %in% examp_beta1, time > 90) %>%
  ggplot(aes(x = time, y = I)) + 
  geom_line() + 
  facet_wrap(vars(beta1), ncol = 1) + 
  theme_bw()
p3 = bind_rows(bifur_reslts_import) %>%
  filter(round(beta1, 5) %in% examp_beta1, time > 94) %>%
  ggplot(aes(x = S, y = I)) + 
  geom_path() + 
  facet_wrap(vars(beta1), ncol = 1) + 
  theme_bw()
cowplot::plot_grid(p1, p2, p3, nrow = 1)

### now run through different vaccination scenarios

vax_import_peak_results = list()
vax_import_trough_results = list()
vax_import_peak_results_grad = list()
vax_import_trough_results_grad = list()
for(i in 1:length(examp_beta1)){
  paras_tmp = paras_importation
  paras_tmp["beta1"] = examp_beta1[i]
  # jump peak
  vax_import_peak_results[[i]] = as.data.frame(
    ode(start, times, seirmod2, paras_tmp, vax_pct = vax_peakstart, rtol = 1e-10, atol = 1e-10)) %>%
    mutate(beta1 = examp_beta1[i], 
           # Re = c(0, diff(Re)), 
           Re = paras_tmp["beta0"] * (1 + paras_tmp["beta1"] * cos(2 * pi * time)) * paras_tmp["sigma"] * S / ((paras_tmp["gamma"] + paras_tmp["mu"]) * (paras_tmp["sigma"] + paras_tmp["mu"])  * paras_tmp["N"]))
  # jump trough
  vax_import_trough_results[[i]] = as.data.frame(
    ode(start, times, seirmod2, paras_tmp, vax_pct = vax_troughstart)) %>%
    mutate(beta1 = examp_beta1[i], 
           Re = paras_tmp["beta0"] * (1 + paras_tmp["beta1"] * cos(2 * pi * time)) * paras_tmp["sigma"] * S / ((paras_tmp["gamma"] + paras_tmp["mu"]) * (paras_tmp["sigma"] + paras_tmp["mu"])  * paras_tmp["N"]))
  # gradual peak
  vax_import_peak_results_grad[[i]] = as.data.frame(
    ode(start, times, seirmod2, paras_tmp, vax_pct = vax_peakstart_grad)) %>%
    mutate(beta1 = examp_beta1[i], 
           Re = paras_tmp["beta0"] * (1 + paras_tmp["beta1"] * cos(2 * pi * time)) * paras_tmp["sigma"] * S / ((paras_tmp["gamma"] + paras_tmp["mu"]) * (paras_tmp["sigma"] + paras_tmp["mu"])  * paras_tmp["N"]))
  # gradual trough
  vax_import_trough_results_grad[[i]] = as.data.frame(
    ode(start, times, seirmod2, paras_tmp, vax_pct = vax_troughstart_grad)) %>%
    mutate(beta1 = examp_beta1[i], 
           Re = paras_tmp["beta0"] * (1 + paras_tmp["beta1"] * cos(2 * pi * time)) * paras_tmp["sigma"] * S / ((paras_tmp["gamma"] + paras_tmp["mu"]) * (paras_tmp["sigma"] + paras_tmp["mu"])  * paras_tmp["N"]))
}

rslts_all_import = bind_rows(vax_import_peak_results_grad) %>%
  mutate(vax_start = "peak", honeymoon_end = "gradual") %>%
  bind_rows(bind_rows(vax_import_trough_results_grad) %>% mutate(vax_start = "trough", honeymoon_end = "gradual")) %>%
  bind_rows(bind_rows(vax_import_trough_results) %>% mutate(vax_start = "trough", honeymoon_end = "jump")) %>%
  bind_rows(bind_rows(vax_import_peak_results) %>% mutate(vax_start = "peak", honeymoon_end = "jump")) %>%
  mutate(stage = ifelse(time < ifelse(vax_start == "peak", time_on_peak, time_on_trough), "pre vax", 
                        ifelse(time < ifelse(vax_start == "peak", time_on_peak, time_on_trough) + nyear_vax, "high vax", "decline vax")), 
         time = ifelse(vax_start == "trough", time - 1, time), 
         # switch peak and trough for beta1 = 0.2 in this case
         vax_start = ifelse(beta1 == 0.2, ifelse(vax_start == "peak", "trough", "peak"), vax_start))


rslts_all_import %>%
  filter(time > 77) %>%
  ggplot() + 
  geom_rect(data = data.frame(t = c(time_on_peak),
                              tend = c(time_on_peak) + nyear_vax,
                              vax_start = c("peak")), 
            aes(xmin = t, xmax = tend, ymin = 0, ymax = Inf), alpha = 0.2) +
  geom_vline(data = data.frame(t = seq(80, 100, 5)), aes(xintercept = t), linetype = "dotted") + 
  geom_line(aes(x = time, y = I, color = stage, linetype = vax_start)) + 
  ggtitle("with importation") + 
  facet_grid(rows = vars(beta1), cols = vars(honeymoon_end)) + 
  scale_color_manual(values = c("blue", "red", "black")) + 
  scale_linetype_manual(values = c("longdash", "solid")) + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank())

rslts_all_import %>%
  filter(time > 77, time < 85) %>%
  ggplot(aes(x = log(S), y = log(I), color = stage, linetype = vax_start)) + 
  geom_path() + 
  ggtitle("with importation") + 
  facet_grid(rows = vars(beta1), cols = vars(honeymoon_end)) + 
  scale_color_manual(values = c("blue", "red", "black")) + 
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank())

rslts_all_import %>%
  filter(time > 77, time < 85) %>%
  ggplot(aes(x = log(Re), y = log(I), color = stage, linetype = vax_start)) + 
  # geom_rect() +
  geom_path() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  ggtitle("with importation") + 
  facet_grid(rows = vars(beta1), cols = vars(honeymoon_end)) + 
  scale_color_manual(values = c("blue", "red", "black")) + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank())

## hmm that changes things, but doesn't seem to have any clearer patterns
## which suggests that it is not small levels of I alone, that are driving the
## transient dynamics