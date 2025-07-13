#### SIR MODEL -----------------------------------------------------------------
#' @param rsv_force list of functions for age-specific rsv_force values
sir_age_structured <- function(t, x, parms, compartments, age_classes, mort, vax_change_times, vax_rates, Fmat,
                                 fert, waifw, adjust_beta_flag = FALSE, print_warnings_flag = FALSE){
  with(as.list(parms),{
    # if(any(abs(t - seq(5, 100, 5)) < 1e-2)){
    #   print(t)
    #   }
    # if(abs(t - 0.5) < 1e-4){browser()}
    # print(t)
    # if(t > 1){browser()}
    x = x[1:(length(x)-4)] # remove beta hat variables for calculations
    if(any(x < 0)){
      if(any(abs(x[which(x<0)]) > 1e-10)){if(print_warnings_flag){print(paste("X NEG! t = ", t)); print(x[which(x < -1e-7)])}}
      x[which(x<0)] = 0}
    nage = length(age_classes)
    ncomp = length(compartments)
    aging <- 1/diff(c(0,age_classes))
    x_mat = matrix(x, ncol = ncomp, nrow = nage, dimnames = list(age_classes, compartments), byrow = TRUE)
    # calculate vaccination rates (assumed at birth)
    if(any(!is.na(vax_change_times))){
      vax_phase = max(which(vax_change_times <= t))
      v = vax_rates[vax_phase]
    }
    else{v = 0}
    # calculate age-specific transmission rates
    beta <- beta0*(beta1*cos(2*pi*(t-p))+1)
    N = sum(x_mat)
    tot_I = sum(x_mat[, "I"])
    tot_S = sum(x_mat[, "S"])
    homogeneous_s = tot_S*diff(c(0, age_classes))/sum(diff(c(0, age_classes))) # account for differences in age bin width
    homogeneous_i = tot_I*diff(c(0, age_classes))/sum(diff(c(0, age_classes)))
    lambda <- waifw%*%(beta/N*(x_mat[, "I"]))
    beta_hat = sum(lambda * x_mat[, "S"])/ ((tot_I * tot_S)/N)
    beta_hat_s = sum(waifw%*%(beta/N*homogeneous_i) * x_mat[, "S"])/ ((tot_I * tot_S)/N) # contribution of heterogeneous susceptibility (i.e., assume constant I)
    beta_hat_i = sum(lambda * homogeneous_s)/ ((tot_I * tot_S)/N) # contribution of heterogeneous infectious (i.e., assume constant S)
    beta_hat_b = sum(waifw%*%(beta/N*homogeneous_i) * homogeneous_s)/ ((tot_I * tot_S)/N) # contribution of heterogeneous mixing (i.e., assume constant S and I)
    # print(beta/beta_hat)
    if(adjust_beta_flag){
      lambda = lambda * beta/beta_hat
    }
    # fertility
    N_fert <- Fmat %*% x
    N_fert <- matrix(N_fert, ncol = ncomp, nrow = nage, dimnames = list(age_classes, compartments))
    N_age <- x_mat * aging
    if(length(age_classes) == 1){
      N_age_in = matrix(0, 1, ncomp, dimnames = list(age_classes, compartments))
      N_age_out = matrix(0, 1, ncomp, dimnames = list(age_classes, compartments))
    }
    else{
      N_age_in <- rbind(rep(0, ncomp), N_age[1:(nrow(N_age)-1), ])
      N_age_out <- rbind(N_age[1:(nrow(N_age)-1),], rep(0, ncomp))
    }
    # calculate age-specific derivatives
    dS <- - (lambda + mort) * x_mat[, "S"] + omega * x_mat[, "R"] +
      N_fert[, "S"] * (1 - v) + N_age_in[, "S"] - N_age_out[, "S"]
    dE <- lambda*x_mat[, "S"] - (mort + sigma) * x_mat[, "E"] +
      N_fert[, "E"] + N_age_in[, "E"] - N_age_out[, "E"]
    dI <- sigma * x_mat[, "E"] - (mort + gamma) * x_mat[, "I"] +
      N_fert[, "I"] + N_age_in[, "I"] - N_age_out[, "I"]
    dR <- gamma * x_mat[, "I"] - (mort + omega) * x_mat[, "R"] +
      N_fert[, "R"] + N_fert[, "S"] * v + N_age_in[, "R"] - N_age_out[, "R"] # NOTE: may want to put v and (1-v in the f matrix instead)
    der <- c(dS, dE, dI, dR)
    names(der) <- sapply(compartments, function(i){paste0(i, "_", age_classes)})
    der <- matrix(c(dS, dE, dI, dR), 
                  nrow = nage, ncol = ncomp)
    der <- c(t(der))
    names(der) <- sapply(age_classes, function(i){paste0(compartments, "_", i)})
    der <- c(der, BH = beta/beta_hat, BHs = beta/beta_hat_s, BHi = beta/beta_hat_i, BHb = beta/beta_hat_b)
    return(list(der))
  })
}

#### HELPERS -------------------------------------------------------------------
#' max_t in weeks
run_ode <- function(age_classes, mort, fert, start_pop, vax_change_times = NA, vax_rates = NA, compartments,
                    waifw = NA, IC_type, IC_manual = NA, max_t, params, beep_flag = FALSE, adjust_beta_flag = FALSE, 
                    print_warnings_flag = FALSE, plot_flag = FALSE, plot_title = NA, dt = 1/12){
  IC = setup_IC(start_pop, age_classes, compartments, mort, fert, IC_type, IC_manual)
  if(any(is.na(waifw))){
    waifw = matrix(1, length(age_classes), length(age_classes)) 
  }
  # run model
  times = seq(0, max_t, dt)
  Fmat <- buildFMatrix(age.classes = age_classes, fert = fert, ncompartments = length(compartments), time.step = 1)
  rslts <- as.data.frame(
    ode(
      y = IC,
      times = times,
      func = sir_age_structured,
      compartments = compartments,
      age_classes = age_classes,
      mort = mort, 
      fert = fert,
      waifw = waifw,
      vax_change_times = vax_change_times, 
      vax_rates = vax_rates,
      parms = params, 
      Fmat = Fmat,
      adjust_beta_flag = adjust_beta_flag, 
      print_warnings_flag = print_warnings_flag
    ))
  if(beep_flag){beep()}
  return(process_results(rslts, plot_flag, plot_title, max_t, dt))
}


buildFMatrix <- function(age.classes=c(1:60, seq(72,120,by=12), seq(180,600,by=60)),  
                         fert =  c(rep(0,66), rep(0.1,7)),
                         ncompartments, time.step, maternal_immunity_flag = FALSE){
  nage <- length(age.classes)
  Fmat <- matrix(0,ncompartments*nage,ncompartments*nage)
  birth_compartment = 1
  for (j in 1:nage) {
    Fmat[birth_compartment,((j-1)*ncompartments+1):(j*ncompartments)] <- rep(fert[j]*time.step, ncompartments)    
  }
  return(Fmat)
}

findStableStruct <- function(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
                             mort=c(rep(1e-9,72),1), 
                             fert =  c(rep(0,66),rep(0.1,7)), time.step = 1){
  nage <- length(age.classes)
  # aging.rate <- time.step/diff(c(0,age.classes))
  aging.rate <- 1/diff(c(0,age.classes))/365 # for daily time steps
  Fmat <- Tmat <- matrix(0,nage,nage)
  for (j in 1:(nage-1)) { 
    Tmat[j,j] <- (1-mort[j]*time.step)*(1-aging.rate[j])
    Tmat[j+1,j] <- (1-mort[j]*time.step)*aging.rate[j]
    Tmat[j,j] <- (1-mort[j])*(1-aging.rate[j])
    Tmat[j+1,j] <- (1-mort[j])*aging.rate[j]
  }
  j <- nage	
  Tmat[j,j] <- (1-mort[j]*time.step)
  Fmat[1,] <- fert*time.step
  # calculate equilibrium values
  stable.age <- Re(eigen(Tmat+Fmat)$vector[,1])
  stable.age <- stable.age/sum(stable.age)
  # adjust for negatives?
  if(any(stable.age < 0)){
    stable.age2 = stable.age
    stable.age[which(stable.age < 0)] = 0
    stable.age[which(stable.age == max(stable.age))] = stable.age[which(stable.age == max(stable.age))] - (sum(stable.age) - 1) 
  }
  lambda <- Re(eigen(Tmat+Fmat)$value[1])
  reprod.value <- Re(eigen(Tmat+Fmat)$vector[1,])
  return(list(stable.age = stable.age, lambda = lambda,
              reprod.value = reprod.value, age.classes = age.classes))
}

# std give 2000 infected individuals across all PA and PB
setup_IC <- function(start_pop, age_classes, compartments, mort, fert, IC_type = "std", IC_manual = NA){
  # indexing - the rows for maternal, susceptible, etc
  indx_comp = rep(compartments, length(age_classes))
  IC <- rep(0,length(age_classes)*length(compartments))
  names(IC) = paste0(rep(compartments, length(age_classes)), "_", 
                     sort(rep(age_classes, length(compartments))))
  if(IC_type == "std"){
    # setup initial conditions
    IC[which(indx_comp == "S")] = start_pop*0.059/length(age_classes) 			 # susceptibles
    IC[which(indx_comp == "I")] = 0.001/length(age_classes)
    IC[which(indx_comp == "R")] = start_pop*0.94/length(age_classes)
  }
  else if(IC_type == "stable-age"){
    expected_stable <- findStableStruct(age.classes = age_classes, mort = mort, fert = fert, time.step = 1)
    IC[which(indx_comp == "S")] = start_pop*0.059*expected_stable$stable.age 			 # susceptibles
    IC[which(indx_comp == "I")] = start_pop*0.001*expected_stable$stable.age 			 # infecteds
    IC[which(indx_comp == "R")] = start_pop*0.94*expected_stable$stable.age 			 # recovereds
    if(any(IC < 0)){browser()}
    if(abs(sum(IC) - start_pop) > 1e-6){browser()}
  }
  # otherwise provide an IC vector (distributed according to stable age distribution)
  else if(IC_type == "manual")
    {
    if(length(IC_manual) == length(compartments)){
      expected_stable <- findStableStruct(age.classes = age_classes, mort = mort, fert = fert, time.step = 1)
      IC[which(indx_comp == "S")] = start_pop*IC_manual["S"]*expected_stable$stable.age
      IC[which(indx_comp == "E")] = start_pop*IC_manual["E"]*expected_stable$stable.age
      IC[which(indx_comp == "I")] = start_pop*IC_manual["I"]*expected_stable$stable.age
      IC[which(indx_comp == "R")] = start_pop*IC_manual["R"]*expected_stable$stable.age
    }
    else if(length(IC_manual == length(IC))){
      IC = IC_manual
    }
  }
  IC = c(IC, BH = 0, BHs = 0, BHi = 0, BHb = 0)
  return(IC)
}

process_results <- function(rslts, plot_flag = FALSE, plot_title = NA, max_t, dt){
  rslts_long <- rslts %>%
    mutate(BH = c(NA, diff(BH))*(1/dt), 
           BHs = c(NA, diff(BHs))*(1/dt), 
           BHi = c(NA, diff(BHi))*(1/dt), 
           BHb = c(NA, diff(BHb))*(1/dt)) %>% # change for different dt
    melt(c("time")) %>%
    tidytable::separate(variable, into = c("variable", "age"), sep = "_") %>%
    mutate(age = as.double(age))
  if(plot_flag){
    rslts_tot <- rslts_long %>% 
      filter(variable == "I") %>%
      summarize(value = sum(value), .by = c("variable", "time"))
    p <- ggplot(data = rslts_tot %>% filter(time > max_t-20), 
                aes(x = time, y = value)) + 
      geom_line() +
      facet_wrap(vars(variable), scales = "free") + 
      theme_bw()
    if(!is.na(plot_title)){
      p <- p + ggtitle(plot_title)
    }
    print(p)
  }
  return(rslts_long)
}

create_polymod_matrix = function(age_classes, plot_flag = FALSE, 
                                 age_classes_to_label = c(seq(12, 60, 12), seq(120, 840, 120))){
  # create polymod using code from Bjornstad book
  data(polymod)
  x = y = polymod$contactor[1:30]
  z = matrix(polymod$contact.rate, ncol = 30, nrow = 30)
  n = length(x)
  # symmetrize
  z2 = (z + t(z))/2
  z3 = as.vector(z2)
  xy = data.frame(x = rep(x[1:n], n), y = rep(y[1:n], each = n))
  polysmooth = Tps(xy, z3, df = 100)
  # surface(polysmooth, xlab = "", ylab = "", col = gray((12:32)/32))
  # annualize & symmetrize
  ps = predict(polysmooth, x = expand.grid(age_classes, age_classes))
  ps2 = matrix(ps, ncol = length(age_classes))
  ps2 = ps2 + t(ps2)
  ps2[ps2<0] = 0 # EH ADDED: IS THIS OKAY?
  W = ps2/mean(ps2)
  if(plot_flag){
    # plot W matrix
    p <- ggplot(data = melt(W), aes(x = Var1, y = Var2, fill = value)) + 
      geom_tile() + 
      scale_fill_viridis_c() + 
      scale_x_continuous(expand = c(0,0),
                         breaks = which(age_classes %in% age_classes_to_label),
                         labels = age_classes_to_label,
                         name = "age (years)") +
      scale_y_continuous(expand = c(0,0),
                         breaks = which(age_classes %in% age_classes_to_label),
                         labels = age_classes_to_label,
                         name = "age (years)")
    print(p)
  }
  return(W)
}

#' should be list
match_on_age_or_inc = function(age_classes, mort, fert, start_pop, compartments, 
                               vax_pct = 0, params, waifws, max_t = 100, plot_flag = FALSE){
  # run ODE with flat WAIFW
  flat = run_ode(age_classes = age_classes, mort = mort, fert = fert, 
                 start_pop = start_pop, compartments = compartments, params = params, 
                 vax_change_times = c(0), vax_rates = c(vax_pct),
                 IC_type = "std", max_t = max_t, dt = 1/12)
  # get age distribution of cases
  age_dist_target = flat %>% 
    filter(time == max(time), variable %in% c("I")) %>% 
    arrange(age) %>%
    pull(value)
  # loop through beta multipliers for all other waifw 
  waifw_scalars = vector("list", length(waifws))
  for(i in 1:length(waifws)){
    print(paste0("starting waifw ", i, "/", length(waifws)))
    waifw_scalars[[i]] = match_one_waifw(age_classes = age_classes, mort = mort, fert = fert, 
                    start_pop = start_pop, compartments = compartments, params = params,
                    vax_pct = vax_pct,max_t = max_t,
                    age_dist_flat = age_dist_target, new_waifw = waifws[[i]])
  }
  waifw_scalars2 = bind_rows(waifw_scalars, .id = "waifw_id") %>%
    mutate(waifw_id = as.integer(waifw_id) + 1) %>% # assuming in the code, waifw1 = flat (not included in this scalaing analysis)
    melt(c("scalar", "waifw_id")) %>%
    mutate(best_value = min(abs(value), na.rm = TRUE), .by = c("waifw_id", "variable")) %>%
    mutate(best_value = ifelse(abs(value) == best_value, TRUE, FALSE))
  if(plot_flag){
    p = ggplot(data = waifw_scalars2, 
           aes(x = scalar, y = value)) + 
      geom_line() + 
      geom_point(aes(color = best_value)) + 
      facet_grid(cols = vars(waifw_id), rows = vars(variable), scales = "free") + 
      ggtilte(paste0("v = ", vax_pct)) +
      scale_color_manual(values = c("black", "red")) + 
      theme_bw()
    print(p)
  }
  return(waifw_scalars2)
}

match_one_waifw = function(age_classes, mort, fert, start_pop, compartments, params, 
                           age_dist_flat, new_waifw, vax_pct, max_t, 
                           scalars = exp(seq(-3, 3, 0.25))){
  bin_width = diff(c(0, age_classes))
  flat_mean_age = sum(age_dist_flat*(age_classes-bin_width/2))/sum(age_dist_flat)
  new_waifw_results = data.frame(scalar = scalars, tot_I_diff = NA, mean_age_diff = NA) #abs_age_diff = NA, 
  for(i in 1:length(scalars)){
    tmp_pars = params
    tmp_pars["beta0"] = tmp_pars["beta0"]*scalars[i]
    tmp_out = tryCatch(
      #this is the chunk of code we want to run
      {run_ode(age_classes = age_classes, mort = mort, fert = fert, 
               start_pop = start_pop, compartments = compartments, params = tmp_pars,
               vax_change_times = c(0), vax_rates = c(vax_pct),
               waifw = new_waifw, IC_type = "std", max_t = max_t, dt = 1/12)
        #when it throws an error, the following block catches the error
      }, error = function(msg){
        return(data.frame(time = NA))
      })
    if(any(is.na(tmp_out$time))){
      print(paste0("error in scalar of ", scalars[i]))
      # new_waifw_results[i, "abs_age_diff"] = NA
      new_waifw_results[i, "tot_I_diff"] = NA
      new_waifw_results[i, "mean_age_diff"] = NA
    }
    else{
      age_dist_new_waifw = tmp_out %>%
        filter(time == max(time), variable == "I") %>% 
        arrange(age) %>%
        pull(value)
      new_waifw_mean_age = sum(age_dist_new_waifw*(age_classes - bin_width/2))/sum(age_dist_new_waifw)
      # new_waifw_results[i, "abs_age_diff"] = sum(abs(age_dist_new_waifw - age_dist_flat))
      new_waifw_results[i, "tot_I_diff"] = sum(age_dist_new_waifw) - sum(age_dist_flat)
      new_waifw_results[i, "mean_age_diff"] = new_waifw_mean_age - flat_mean_age
    }
  }
  return(new_waifw_results)
}


