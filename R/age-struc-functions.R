#' @param rsv_force list of functions for age-specific rsv_force values
sir_age_structured <- function(t, x, parms, compartments, age_classes, mort, vax_change_times, vax_rates,
                                 fert, waifw, adjust_beta_flag = FALSE, print_warnings_flag = FALSE){
  with(as.list(parms),{
    x = x[1:(length(x)-3)] # remove beta hat variables for calculations
    # if(any(x < 0)){print("X NEG!")}
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
    # print(beta/beta_hat)
    if(adjust_beta_flag){
      lambda = lambda * beta/beta_hat
    }
    # fertility
    Fmat <- buildFMatrix(age.classes = age_classes, fert = fert, ncompartments = ncomp, time.step = 1)
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
    # browser()
    der <- c(der, BH = beta/beta_hat, BHs = beta/beta_hat_s, BHi = beta/beta_hat_i)
    return(list(der))
  })
}

vax_change_times#' max_t in weeks
run_ode <- function(age_classes, mort, fert, start_pop, vax_change_times = NA, vax_rates = NA, 
                    waifw = NA, IC_type, max_t, params, beep_flag = FALSE, adjust_beta_flag = FALSE, 
                    print_warnings_flag = FALSE, plot_flag = FALSE, plot_title = NA){
  IC = setup_IC(start_pop, age_classes, compartments, mort, fert, IC_type)
  if(any(is.na(waifw))){
    waifw = matrix(1, length(age_classes), length(age_classes)) 
  }
  # run model
  times = seq(0, max_t, 1/365)
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
      adjust_beta_flag = adjust_beta_flag, 
      print_warnings_flag = print_warnings_flag
    ))
  if(beep_flag){beep()}
  return(process_results(rslts, plot_flag, plot_title, max_t))
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
  aging.rate <- 1/diff(c(0,age_classes))/365 # for daily time steps
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
  lambda <- Re(eigen(Tmat+Fmat)$value[1])
  reprod.value <- Re(eigen(Tmat+Fmat)$vector[1,])
  return(list(stable.age = stable.age, lambda = lambda,
              reprod.value = reprod.value, age.classes = age.classes))
}

# std give 2000 infected individuals across all PA and PB
setup_IC <- function(start_pop, age_classes, compartments, mort, fert, type = "std"){
  # indexing - the rows for maternal, susceptible, etc
  indx_comp = rep(compartments, length(age_classes))
  if(type == "std"){
    # setup initial conditions
    IC <- rep(0,length(age_classes)*length(compartments))
    names(IC) = paste0(rep(compartments, length(age_classes)), "_", 
                       sort(rep(age_classes, length(compartments))))
    IC[which(indx_comp == "S")] = start_pop*0.059/length(age_classes) 			 # susceptibles
    IC[which(indx_comp == "I")] = 0.001/length(age_classes)
    IC[which(indx_comp == "R")] = start_pop*0.94/length(age_classes)
    IC = c(IC, BH = 0, BHs = 0, BHi = 0)
    return(IC)
  }
  if(type == "stable-age"){
    IC <- rep(0,length(age_classes)*length(compartments))
    names(IC) = paste0(rep(compartments, length(age_classes)), "_", 
                       sort(rep(age_classes, length(compartments))))
    expected_stable <- findStableStruct(age.classes = age_classes, mort = mort, fert = fert, time.step = 1)
    IC[which(indx_comp == "S")] = start_pop*0.059*expected_stable$stable.age 			 # susceptibles
    IC[which(indx_comp == "I")] = start_pop*0.001*expected_stable$stable.age 			 # infecteds
    IC[which(indx_comp == "R")] = start_pop*0.94*expected_stable$stable.age 			 # recovereds
    if(any(IC < 0)){browser()}
    if(abs(sum(IC) - start_pop) > 1e-6){browser()}
    IC = c(IC, BH = 0, BHs = 0, BHi = 0)
    return(IC)
  }
}

process_results <- function(rslts, plot_flag = FALSE, plot_title = NA, max_t){
  rslts_long <- rslts %>%
    mutate(BH = c(NA, diff(BH))*365, 
           BHs = c(NA, diff(BHs))*365, 
           BHi = c(NA, diff(BHi))*365) %>% # change for different dt
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



