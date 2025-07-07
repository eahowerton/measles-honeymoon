library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidytable)
library(epimdr2)
library(cowplot)
library(dplyr)
library(beepr)

source("R/age-struc-functions.R")

#### MODEL SETUP ---------------------------------------------------------------
compartments = c("S", "E", "I", "R")
times = seq(0, 100, by = 1/365)
times_long = seq(0, 300, by = 1/365)

no_vax = approxfun(times_long, rep(0, length(times_long)))

paras = c(mu = 1/50, N = 1, beta0 = 1000, beta1 = 0.2, omega = 0,
          sigma = 365/8, gamma = 365/5, vax_pct = 0, delta = 0, p = 0)

start_pop = 1
fert = paras["mu"]
mort = paras["mu"]

#### TEST WITH ONE COMPARTMENT -------------------------------------------------
age_classes1 =  c(100)
rslts1 <- run_ode(
  age_classes = age_classes1, mort = mort, fert = fert, start_pop = start_pop, 
  IC_type = "std", max_t = max(times_long), params = paras,
  plot_flag = TRUE, plot_title = "1 age class")

#### TEST WITH TWO COMPARTMENTS ------------------------------------------------
age_classes2 = c(20, 40)
rslts2 <- run_ode(
  age_classes = age_classes2, mort = rep(mort, length(age_classes2)),
  fert = rep(fert, length(age_classes2)), start_pop = start_pop,
  IC_type = "std", max_t = max(times_long), params = paras,
  plot_flag = TRUE, plot_title = "2 age classes")

#### TEST WITH REALISTIC STRUCTURE ---------------------------------------------
# up to 5 years in months, up to 10 in years, and up to 70 in 5 years 
# (multiply by 4 to go from months to weeks)
age_classes = c(seq(1, 5, 1/2), 6:9, seq(10,70,by=10))
bin_width = diff(c(0, age_classes))
names(bin_width) = age_classes

start.time <- Sys.time()
rslts3 <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop,
  IC_type = "std", max_t = max(times_long), params = paras, beep_flag = TRUE,
  plot_flag = TRUE, plot_title = "realistic structure, constant WAIFW")
Sys.time() - start.time

#### NOW ADD VACCINATION -------------------------------------------------------
vax_change_times = c(0, 190, 195)
vax_rates = c(0, 0.8, 0.4)

start.time <- Sys.time()
rslts3_v <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop, 
  vax_change_times = vax_change_times, vax_rates = vax_rates,
  IC_type = "std", max_t = max(times_long), params = paras, beep_flag = TRUE,
  plot_flag = TRUE, plot_title = "realistic structure, constant WAIFW")
Sys.time() - start.time

#### TEST WITH POLYMOD CONTACTS ------------------------------------------------
# now add more realistic mixing with POLYMOD data
W = create_polymod_matrix(age_classes)

start.time <- Sys.time()
rslts4_u <- run_ode(
  age_classes = age_classes, mort = rep(mort, length(age_classes)),
  fert = rep(fert, length(age_classes)), start_pop = start_pop, waifw = W,
  IC_type = "std", max_t = max(times_long)+75, params = paras, beep_flag = TRUE,
  adjust_beta_flag = TRUE, plot_flag = TRUE, plot_title = "realistic structure, POLYMOD (adjust beta)")
Sys.time() - start.time


