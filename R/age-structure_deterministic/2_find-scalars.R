library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidytable)
library(epimdr2)
library(cowplot)
library(dplyr)
library(beepr)
library(fields)
library(rootSolve)

source("R/age-structure_deterministic/0_helper-functions.R")
source("R/age-structure_deterministic/0_SIR-age.R")

#### MODEL SETUP ---------------------------------------------------------------
source("R/age-structure_deterministic/1_setup-WAIFW.R") # this will add waifw to environment, which is list of WAIFW matrices to test
n_waifw = length(waifw)

compartments = c("S", "I", "R")

R0 <- 14
paras = c(mu = 1/50, N = 500000, beta1 = 0,
          gamma = 365/14, delta = 0, p = 0) # delta = 1e-4
paras["beta0"]  = (paras["gamma"] + paras["mu"])*R0

# R0
with(as.list(paras), beta0/(gamma + mu))

start_vax = c(0.9, 0.92, 0.94, 0.96)

#### GET R0 VALUES FOR EACH WAIFW ----------------------------------------------
stable_age = findStableStruct(age_classes, mort, fert)$stable.age

get_Rt(waifw[[1]], stable_age, paras["beta0"], paras["gamma"], paras["N"])

#### FIND SCALARS --------------------------------------------------------------
# each WAIFW matrix yields a different equilibrium number of infections 
scalars_v = vector("list", length(start_vax))
for(i in 1:length(start_vax)){
  print(paste0("start vax = ", start_vax[i]))
  scalars_v[[i]] = match_on_age_or_inc(age_classes = age_classes, mort = mort, 
                                       fert = fert, start_pop = paras["N"], 
                                       compartments = compartments, 
                                       vax_pct = start_vax[i], max_t = 200, 
                                       params = paras, waifws = waifw[-1])
}

scalars_v_full =  bind_rows(scalars_v, .id = "vax_id") %>%
  mutate(vax_id = as.integer(vax_id)) %>%
  left_join(data.frame(vax_id = 1:length(start_vax), 
                       v = start_vax))

saveRDS(scalars_v_full, "data/output-data/scalars_by_v.rds")
