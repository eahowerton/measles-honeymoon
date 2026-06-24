#### ANALYSIS OF REDUCTIONS IN VAX COVERAGE IN SIMPLE SIR MODEL

# First, we make some simplifying assumptions
# 1. assume that the population has been under a constant vax program for 
#    sufficiently long that the population is at a vax equilibrium
# 2. assume that measles is locally extinct
# 
# then, the question is when should we expect a rebound outbreak
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

# honeymoon time derived analytically
return_time = equilibrium_S %>%
  mutate(honeymoon_time = -1/mu * log(((1/R0) - (1-new_v))/(S_star - (1-new_v))))

p = return_time %>%
  ggplot(aes(x = new_v, y = honeymoon_time, color = as.factor(v))) + 
  geom_path(linewidth = 1) +
  scale_color_manual(values = RColorBrewer::brewer.pal(5, "Greys")[2:4], name = "pre-decline\nimmunization", 
                     labels = paste0(pre_v*100, "%")) +
  scale_x_continuous(breaks = seq(0, 9, 0.3), expand = c(0,0),
                     label = scales::percent, name = "post-decline immunization") + 
  scale_y_continuous(name = "theoretical honeymoon time\n(time to Re > 1)") + 
  theme_bw(base_size = 7) + 
  theme(legend.position = "bottom")
ggsave("figures/honeymoon_SIR_full.pdf", p, width = 4, height = 4.5)

# test how this is sensitive to birth rate
broader_return_time = expand.grid(mu = c(1/50, 1/80), 
                            beta0 = paras["beta0"], 
                            v = seq(0.93, 1, 0.005), 
                            gamma = paras["gamma"], 
                            new_v = seq(0, 0.90, 0.01)) %>%
  mutate(R0 = (beta0)/((gamma + mu)), 
         herd_imm = 1-1/R0,
         S_star = ifelse(v > herd_imm, 1-v,
                         (mu + gamma)/beta0)) %>%
  mutate(honeymoon_time = -1/mu * log(((1/R0) - (1-new_v))/(S_star - (1-new_v))))

# overlay manual calculations to double check math
# contours_manual = expand.grid(
#   new_v = seq(0, 0.9, 0.01),
#   mu = 1/50,
#   R0 = 17,
#   honeymoon_time = c(1,5,10,15)
# ) %>%
#   mutate(S_star = 1/R0,
#          line = exp(honeymoon_time*mu)*(1 - S_star) + (1-exp(mu*honeymoon_time))*new_v)

mu_labs = paste0("'birth rate' == 1/", c(80, 50), "~(years^-1)")
names(mu_labs) = c(-80,-50)

g = ggplot(broader_return_time %>% mutate(mu_lab = -1/mu), aes(x = new_v, y = v)) + 
  geom_tile(aes(fill = honeymoon_time)) + 
  geom_contour(aes(z = honeymoon_time), color = "black", breaks = c(1, 5, 10, 15), linewidth = 0.25) + 
  metR::geom_text_contour(aes(z = honeymoon_time), stroke = 0.15, breaks = c(1, 5, 10, 15), size = 2) +
  # geom_line(data = contours_manual, aes(x = new_v, y = line, group = honeymoon_time), linetype = 'dashed', color = "red") +
  facet_wrap(vars(mu_lab), ncol = 1, labeller = as_labeller(mu_labs, default = label_parsed)) +#, switch = "y") + 
  scale_fill_distiller(palette = "Blues", direction = 0, name = "theoretical\nhoneymoon\ntime") +
  scale_linetype_manual(values = c("longdash", "solid", "dotted")) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0, 0.9, 0.3), label = scales::percent, 
                     name = "post-decline immunization") + 
  scale_y_continuous(expand = c(0,0), label = scales::percent, 
                     name = "pre-decline immunization") + 
  theme_bw(base_size = 7) + 
  theme(legend.position = "bottom", 
        strip.background = element_blank())#, 
        # strip.placement = "outside")

cowplot::plot_grid(p + theme(legend.position = c(0.15, 0.77)), 
                   g + theme(legend.position = "right"), 
                   #align = "h", axis = "tb", 
                   rel_widths = c(0.55, 0.45),
                   labels = c("A", "B"), label_size = 8)
ggsave("figures/honeymoon_SIR_with_linear_heatmap.pdf", width = 6, height = 3)

### values for text
return_time %>%
  filter(v == 0.95, new_v %in% c(0.8, 0.75, 0.25, 0.3))

return_time %>%
  filter(new_v == 0)

