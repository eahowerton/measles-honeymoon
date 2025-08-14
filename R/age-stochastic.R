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

source("R/age-struc-functions.R")
# okay let's try a stochastic model to see how frequently importations into the unvax group occur?

# define two patches of individuals, patch 1 and patch 2
# we conceptualize (losely) of patch 1 as the general population, and patch 2 as
# the vaccine hesistant population

# then, we can define a few values
# Ni = number of inidividuals in patch i
# vi = vax coverage in patch i
# v = sum(vi*Ni)/sum(Ni) = (v1*N1 + v2*N2)/(N1 + N2)
# part of the challenge is you can get the same v by varying v1/v2 and N1/N2
# hmm, it would be nice to have a single metric (from 0 to 1) that represents
# how separated/connected these two groups are

# but wait, can we simplify this to just an initial condition? 
# so let's say that there are these two groups
# and there is an initial infection into the general population
# then, the force of infection in the smaller population will be 
# lambda_2 = (beta*pi*I1 + beta*I2)*S2/N2 where I1 = 1

# okay, idea for tomorrow: 
# 1. stochastic simulation with random importation time (or 100 runs of importations at time t, scanning across a range of t)
# 2. calculate outbreak size in normal SIR model
# 3. then, somehow add some differential mixing between them and see how it changes
# 4. and repeat with age structure? (maybe try the deterministic without age too...)



#' @param x array of state values (compartments x sims)
#' @param delta_t double, length of time step
#' @param n_sims integer, number of stochastic replicates
#' @param compartments vector of model compartments names with age breakdowns
#' @param age_classes vector of upper limits on age classes
#' @param params vector of parameter values
#' @param waifw contact matrix
#' @param import_times vector of times, if t %in% import_times, introduce an infection
stoch_step <- function (x, t, delta_t, n_sims, compartments, age_classes, params, waifw, Fmat, import_times) {
  with(as.list(params), {
    ## some setup
    n_compartments = length(compartments)
    n_age = length(age_classes)
    n_reps = n_sims * n_age
    S_indx = which(substr(compartments,1,1) == "S")
    I_indx = which(substr(compartments,1,1) == "I")
    R_indx = which(substr(compartments,1,1) == "R")
    # demography transitions (NOTE: THESE COME FIRST...)
    aging <- 1/diff(c(0,age_classes))
    aging <- aging[sort(rep(1:length(aging), n_compartments/n_age))] # repeat for each compartment and order correctly
    N_age <- rbinom(n = prod(dim(x)), size = x, prob = 1-exp(-aging*delta_t))
    N_age <- matrix(N_age, nrow = n_compartments, ncol = n_sims)
    N_age_in <- rbind(matrix(rep(0, n_sims*3), ncol = n_sims), 
                      N_age[-((nrow(N_age)-2):nrow(N_age)), ])
    N_age_out <- N_age # all death = age out of final compartment
    # fertility (by vax type)
    # assume (for now) that all who age out of final compartment are born again
    N_fert <- colSums(N_age[((nrow(N_age)-2):nrow(N_age)), ])
    N_fert_vax = matrix(0, nrow = n_compartments, ncol = n_sims)
    N_fert_vax[3, ] = rbinom(n = length(N_fert), size = N_fert, prob = v) # vaccinated go into R_0
    N_fert_vax[1, ] = N_fert - N_fert_vax[3, ]                            # the rest (not vaccinated) go into S_0 
    x_demog <- x + N_fert_vax + N_age_in - N_age_out
    # get force of infection
    beta <- beta0*(beta1*cos(2*pi*(t-p))+1)
    N <- colSums(x)
    lambda <- waifw%*%(beta/N*x_demog[I_indx, ])
    StoI = rbinom(n = n_reps, size = x_demog[S_indx, ], prob = 1-exp(-lambda*delta_t))
    ItoR = rbinom(n = n_reps, size = x_demog[I_indx, ], prob = 1-exp(-gamma*delta_t))
    imports = rbinom(n = n_reps, size = x_demog[S_indx, ], prob = 1-exp(-delta*delta_t))
    # imports = as.matrix(imports, nrow = n_compartments, ncol = n_sims)
    if(round(t,4) %in% import_times){ # import in age 5 (for now)
      age_5_indx = which(age_classes == 5) + length(age_classes)*(0:(n_sims-1))
      imports[age_5_indx] = imports[age_5_indx] + 1
    }
    # implement transitions
    new_x = x_demog
    new_x[S_indx, ] <- new_x[S_indx, ] - StoI - imports
    new_x[I_indx, ] <- new_x[I_indx, ] + StoI + imports - ItoR
    new_x[R_indx, ] <- new_x[R_indx, ] + ItoR
    return(new_x)
  })
}

# then we need to wrap this in a function that iterates across time
run_sims = function(IC, params, n_timesteps, delta_t, compartments, age_classes, waifw, import_times){
  with(as.list(params),{
    Fmat <- buildFMatrix_withv(age.classes = age_classes, fert = fert, 
                         ncompartments = length(unique(substr(compartments, 1, 1))), 
                         time.step = delta_t, v = v)
    n_sims = ncol(IC)
    ret = array(integer(), c(length(compartments), n_sims, n_timesteps), 
                dimnames = list(NULL, names(IC), NULL))
    # ret dimensions - 1: sims, 2: classes, 3: time
    ret[,,1] = IC
    # run first phase of intervention (distanced, no testing)
    t = 0
    for(ts in 2:n_timesteps){
      t = t + delta_t
      ret[,,ts] <- stoch_step(x = ret[,, ts - 1], t = t, delta_t = delta_t,
                              n_sims = n_sims, compartments = compartments,
                              age_classes = age_classes, params = params, waifw = waifw, Fmat = Fmat, import_times)
    }
    return(ret)
  })
}
  
buildFMatrix_withv <- function(age.classes=c(1:60, seq(72,120,by=12), seq(180,600,by=60)),  
                         fert =  c(rep(0,66), rep(0.1,7)), v,
                         ncompartments, time.step, maternal_immunity_flag = FALSE){
  nage <- length(age.classes)
  Fmat <- matrix(0,ncompartments*nage,ncompartments*nage)
  unvax_compartment = 1
  vax_compartment = 3
  for (j in 1:nage) {
    Fmat[unvax_compartment,((j-1)*ncompartments+1):(j*ncompartments)] <- rep(fert[j]*time.step*(1-v), ncompartments)
    Fmat[vax_compartment,((j-1)*ncompartments+1):(j*ncompartments)] <- rep(fert[j]*time.step*v, ncompartments)    
  }
  return(Fmat)
}

# start with deterministic equilibrium of pre-vax values
# just do stochastic rebound
source("R/1_setup-WAIFW.R") # this will add waifw to environment, which is list of WAIFW matrices to test
n_waifw = length(waifw)

paras = c(mu = 1/50, N = 500000, beta0 = 365, beta1 = 0, 
          gamma = 365/14, v = 0.2, p = 0, delta = 0)

n_sims = 1000
dt = 1/52
start_vax = 0.95
rebound_vax = 0.6
paras_rebound = paras
paras_rebound["v"] = rebound_vax

IC_manual = c(1-start_vax, 0, start_vax)
names(IC_manual) = c("S", "I", "R")
IC = setup_IC(start_pop = paras["N"], age_classes, c("S", "I", "R"), fert = fert, 
              mort = mort, IC_type = "manual", IC_manual = IC_manual, dt)
IC = round(IC[1:(length(IC)-5)])
IC[which(IC == max(IC))] = IC[which(IC == max(IC))] - 100
IC[length(IC)] = 101
IC_mat = matrix(rep(IC, n_sims), ncol = n_sims)

outbreak_thresh = 50 
tst_import_times = seq(0.5, 10.5, 0.5)
sim_import_times = vector("list", length(tst_import_times))
for(i in 1:length(tst_import_times)){
  # print(paste0("import time ", i, "/", length(tst_import_times)))
  tmp = run_sims(IC_mat, paras_rebound, n_timesteps = (2 + tst_import_times[i])/dt, delta_t = dt, 
                                   compartments = names(IC), age_classes, waifw = waifw[[1]], import_time = tst_import_times[i]) %>%
    melt() %>%
    rename(compartment_id = Var1, 
           sim_id = Var2, 
           time = Var3) %>%
    left_join(data.frame(compartment_id = 1:length(names(IC)), 
                         compartment = names(IC))) %>%
    mutate(variable = substr(compartment, 1,1), 
           age = as.double(substr(compartment, 3, length(compartment))), 
           time = time*dt) 
  tmp_tot = tmp %>% 
    summarize(value = sum(value), .by = c("variable", "sim_id", "time")) 
  print(
    ggplot(data = tmp_tot %>% filter(variable == "I"), # just plot first 50 to get an idea
           aes(x = time, y = value, group = sim_id)) + 
      geom_line(alpha = 0.2) + 
      ggtitle(paste0("import time = ", tst_import_times[i])) +
      theme_bw()
  )
  pct_outbreaks = tmp_tot %>%
    # now calculate outbreak pct and only save this since doing more runs
    filter(variable == "I") %>%
    summarize(outbreak_flag = ifelse(any(value > outbreak_thresh), 1, 0), .by = c("sim_id")) %>%
    summarize(pct_outbreaks = sum(outbreak_flag)/n())
  tot_s = tmp_tot %>% filter(variable == "S", round(time,4) == tst_import_times[i])
  sim_import_times[[i]] = data.frame(import_time = tst_import_times[i], 
                                     pct_outbreaks = pct_outbreaks,
                                     quantile = paste0("Q", c(0.05, 0.25, 0.5, 0.75, 0.95)*100),
                                     tot_S = quantile(tot_s$value, c(0.05, 0.25, 0.5, 0.75, 0.95)))
}

p1 = ggplot(data = bind_rows(sim_import_times) %>% select(import_time, pct_outbreaks) %>% unique(),
       aes(x = import_time, y = pct_outbreaks)) + 
  geom_point() + 
  geom_line() + 
  theme_bw()
p2 = bind_rows(sim_import_times) %>%
  mutate(quantile = paste0("Q", quantile*100)) %>%
  dcast(import_time ~ quantile, value.var = "tot_S") %>%
  ggplot(aes(x = import_time)) + 
  geom_segment(aes(xend = import_time, y = Q5, yend = Q95), linewidth = 0.5) +
  geom_segment(aes(xend = import_time, y = Q25, yend = Q75), linewidth = 0.7) +
  geom_point(aes(y = Q50)) +
  geom_line(aes(y = Q50), alpha = 0.5) + 
  theme_bw()
p1/p2


# next: add different release values, different WAIFWs

# next: add assortativity between vax and nonvax
