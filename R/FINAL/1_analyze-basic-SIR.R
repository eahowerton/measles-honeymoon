#### ANALYSIS OF REDUCTIONS IN VAX COVERAGE IN SIMPLE SIR MODEL

# First, we make some simplifying assumptions
# 1. assume that the population has been under a constant vax program for 
#    sufficiently long that the population is at a vax equilibrium
# 2. assume that measles is locally extinct
# 
# then, the question iswhen should we expect a rebound outbreak
# 
# the hypothesis is that this will depend on susceptibility and importation
# susceptibility will depend on birth rates and changes in vaccination
# 
# framing the problem in this way has a few advantages: 
# 1. if we focus on a single importation at different periods of time, we can 
# remove the problem with small numbers of infected individuals.
# in other words, we don't want to see long honeymoons because there I is on the 
# order of e-15...
# 2. different shapes of vaccination curves can just be translated into how they
# affect the build up of susceptible individuals
# 
# so, to test this, we define the "honeymoon period" as time until Re > 1
# this is analytically tractable
# 
# we start without seasonality, and can add it if of interest

#### FUNCTIONS -----------------------------------------------------------------
sirmod = function(t, y, parameters, vax_pct) {
  S = y[1]
  I = y[2]
  R = y[3]
  parameters = c(parameters, v = vax_pct(t))
  with(as.list(parameters), {
    beta = beta0 * (1 + beta1 * cos(2 * pi * (t + p)))
    dS = mu * (1 - v) * N - beta * S * I/N - mu * S - delta * S # additional importations
    dI = beta * S * I/N - (mu + gamma) * I  + delta * S
    dR = mu * N * v + gamma * I - mu * R
    res = c(dS, dI, dR)
    list(res)
  })
}

#### SETUP ---------------------------------------------------------------------
library(dplyr)
library(deSolve)
library(ggplot2)
library(cowplot)
library(reshape2)

# setup parameters for simulation
paras = c(mu = 1/80, N = 1, beta1 = 0, gamma = 365/14, delta = 0, R0 = 17)
paras["beta0"] = paras["R0"] * (paras["gamma"] + paras["mu"])
start = c(S = 0.06, I = 0.001, R = 0.939)
pre_v = c(0.94, 0.95, 0.96)

#### CALCULATE HONEYMOON PERIODS -----------------------------------------------
# first calculate equilibrium susceptible population 
equilibrium_S = expand.grid(mu = paras["mu"], 
                           beta0 = paras["beta0"], 
                           v = pre_v, 
                           gamma = paras["gamma"], 
                           new_v = seq(0, 0.90, 0.01)) %>%
  mutate(R0 = (beta0)/((gamma + mu)), 
         herd_imm = 1-1/R0,
         S_star = ifelse(v > herd_imm, 1-v,
                         (mu + gamma)/beta0))

# so, there will be an outbreak when there is an introduction and 
# Re = S*R0 > 1 , or S* > 1/R0 or S* - 1/R0 > 0
# this quantity (S* - 1/R0) tells us how much "wiggle room" we have, and the
# amount of time until we cross this threshold will depend on the birth rate, and
# the new vaccination rate

# as an example, let's assume that vaccination drops immediately to v_r = 40%
# and stays there. Then, we can calculate the time to rebound based on this new
# vaccination rate where dS/dt = mu * (1 - v_r) - mu * S
# we can solve this for t to get the time until S crosses the threshold making Re > 1

return_time = equilibrium_S %>%
  mutate(time_to_cross = -1/mu * log(((1/R0) - (1-new_v))/(S_star - (1-new_v))))

p_inset = return_time %>%
  filter(new_v < 0.6) %>%
  ggplot(aes(x = new_v, y = time_to_cross, color = as.factor(v))) + 
  geom_path(size = 1) +
  scale_color_manual(values = RColorBrewer::brewer.pal(5, "Greys")[2:4], name = "coverage\nbefore drop", 
                     labels =  paste0(pre_v*100, "%")) +
  scale_x_continuous(label = scales::percent, name = "coverage after drop") +
  scale_y_continuous(name = "time to Re > 1") + 
  theme_bw() + 
  theme(legend.position = "bottom")
ggsave("R/FINAL/figures/honeymoon_SIR_sub.pdf", p_inset, width = 4, height = 4.5)

p = return_time %>%
  ggplot(aes(x = new_v, y = time_to_cross, color = as.factor(v))) + 
  geom_path(size = 1) +
  scale_color_manual(values = RColorBrewer::brewer.pal(5, "Greys")[2:4], name = "coverage\nbefore drop", 
                     labels = paste0(pre_v*100, "%")) +
  scale_x_continuous(breaks = seq(0, 9, 0.3), expand = c(0,0),
                     label = scales::percent, name = "coverage after drop") + 
  scale_y_continuous(name = "time to Re > 1") + 
  theme_bw(base_size = 7) + 
  theme(legend.position = "bottom")


# ggdraw(p + theme_half_open(12)) +
#   draw_plot(inset, .45, .45, .5, .5) +
#   draw_plot_label(
#     c("A", "B"),
#     c(0, 0.45),
#     c(1, 0.95),
#     size = 12
#   )
# 
ggsave("R/FINAL/figures/honeymoon_SIR_full.pdf", width = 4, height = 4.5)

# so there is non-linearity in how long we expect the honeymoon to be, which 
# become more extreme as vaccination coverage gets higher
# in other words, preventing a 10% drop in vax coverage buys you more time if it's
# not dropping too low, and this effect is stronger if your vaccination coverage
# was high to start with 


broader_return_time = expand.grid(mu = c(1/50, 1/80), 
                            beta0 = paras["beta0"], 
                            v = seq(0.93, 1, 0.005), 
                            gamma = paras["gamma"], 
                            new_v = seq(0, 0.90, 0.01)) %>%
  mutate(R0 = (beta0)/((gamma + mu)), 
         herd_imm = 1-1/R0,
         S_star = ifelse(v > herd_imm, 1-v,
                         (mu + gamma)/beta0)) %>%
  mutate(time_to_cross = -1/mu * log(((1/R0) - (1-new_v))/(S_star - (1-new_v))))

g = ggplot(broader_return_time, aes(x = new_v, y = v)) + 
  geom_tile(aes(fill = time_to_cross)) + 
  geom_contour(aes(z = time_to_cross), color = "black", breaks = c(1, 5, 10, 15), linewidth = 0.25) + 
  metR::geom_text_contour(aes(z = time_to_cross), stroke = 0.15, breaks = c(1, 5, 10, 15), size = 2) +
  facet_wrap(vars(-1/mu), ncol = 1, labeller = label_both) +#, switch = "y") + 
  scale_fill_distiller(palette = "Blues", direction = 0, name = "time to\nRe > 1") +
  scale_linetype_manual(values = c("longdash", "solid", "dotted")) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0, 0.9, 0.3), label = scales::percent, 
                     name = "coverage after drop") + 
  scale_y_continuous(expand = c(0,0), label = scales::percent, 
                     name = "coverage before drop") + 
  theme_bw(base_size = 7) + 
  theme(legend.position = "bottom", 
        strip.background = element_blank())#, 
        # strip.placement = "outside")

cowplot::plot_grid(p + theme(legend.position = c(0.15, 0.77)), 
                   g + theme(legend.position = "right"), 
                   #align = "h", axis = "tb", 
                   rel_widths = c(0.55, 0.45),
                   labels = c("A", "B"), label_size = 8)
ggsave("R/FINAL/figures/Figure2.pdf", width = 6, height = 3)

