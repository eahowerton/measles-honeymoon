#### SIR MODEL -----------------------------------------------------------------
sir_age_structured = function(t, x, parms, compartments, age_classes, mort, vax_change_times, vax_rates, Fmat,
                               fert, waifw, adjust_beta_flag = FALSE, print_warnings_flag = FALSE, burnin = 0){
  with(as.list(parms),{
    x = x[1:(length(x)-4)] # remove beta hat variables for calculations
    if(any(x < 0)){x[which(x<0)] = 0; if(t< burnin){x[seq(2, length(x),3)] = x[seq(2, length(x),3)] + 10}; print(paste0("adj at t = ", t))} # check x is non-neg
    nage = length(age_classes)
    ncomp = length(compartments)
    aging = 1/diff(c(0,age_classes))
    x_mat = matrix(x, ncol = ncomp, nrow = nage, dimnames = list(age_classes, compartments), byrow = TRUE)
    # calculate vaccination rates (assumed at birth)
    if(any(!is.na(vax_change_times))){
      vax_phase = max(which(vax_change_times <= t))
      v = vax_rates[vax_phase]
    }
    else{v = 0}
    # calculate age-specific transmission rates
    beta = beta0*(beta1*cos(2*pi*(t-p))+1)
    N = sum(x_mat)
    tot_I = sum(x_mat[, "I"])
    tot_S = sum(x_mat[, "S"])
    # account for differences in age bin width
    homogeneous_s = tot_S*diff(c(0, age_classes))/sum(diff(c(0, age_classes))) 
    homogeneous_i = tot_I*diff(c(0, age_classes))/sum(diff(c(0, age_classes)))
    lambda = waifw %*% (beta/N * x_mat[, "I"])
    beta_hat = sum(lambda * x_mat[, "S"])/ ((tot_I * tot_S)/N)
    if(adjust_beta_flag){
      lambda = lambda * beta/beta_hat
    }
    # fertility
    N_fert = Fmat %*% x
    N_fert = matrix(N_fert, ncol = ncomp, nrow = nage, dimnames = list(age_classes, compartments))
    N_age = x_mat * aging
    # aging
    if(length(age_classes) == 1){
      N_age_in = matrix(0, 1, ncomp, dimnames = list(age_classes, compartments))
      N_age_out = matrix(0, 1, ncomp, dimnames = list(age_classes, compartments))
    }
    else{
      N_age_in = rbind(rep(0, ncomp), N_age[1:(nrow(N_age)-1), ])
      N_age_out = rbind(N_age[1:(nrow(N_age)-1),], rep(0, ncomp))
    }
    # calculate age-specific derivatives
    dS = - lambda * x_mat[, "S"] - delta * x_mat[, "S"] - mort * x_mat[, "S"] +
      N_fert[, "S"] * (1 - v) + N_age_in[, "S"] - N_age_out[, "S"]
    dI = lambda * x_mat[, "S"] + delta * x_mat[, "S"] - gamma * x_mat[, "I"] - mort * x_mat[, "I"] +
      N_fert[, "I"] + N_age_in[, "I"] - N_age_out[, "I"]
    dR = gamma * x_mat[, "I"] - mort * x_mat[, "R"] +
      N_fert[, "R"] + N_fert[, "S"] * v + N_age_in[, "R"] - N_age_out[, "R"] 
    # NOTE: may want to put v and (1-v in the f matrix instead)
    der = c(dS, dI, dR)
    names(der) = sapply(compartments, function(i){paste0(i, "_", age_classes)})
    der = matrix(c(dS, dI, dR), 
                  nrow = nage, ncol = ncomp)
    der = c(t(der))
    names(der) = sapply(age_classes, function(i){paste0(compartments, "_", i)})
    # browser()
    der = c(der, 
             BH = sum(lambda*x_mat[, "S"]),
             BHs = sum(lambda*homogeneous_s), 
             BHi = sum(waifw%*%(beta/N*homogeneous_i)*x_mat[, "S"]), 
             BHb = sum(waifw%*%(beta/N*homogeneous_i)*homogeneous_s))
    return(list(der))
  })
}


