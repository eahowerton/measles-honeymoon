#### AGE-STRUCTURED SIR MODEL ---------------------------------------------------
#' age structured SIR model (to pass to deSolve::ode())
#' @param t time
#' @param x current states (all ages and compartments)
#' #' @param params named vector of parameter values
#' @param compartments vector of ODE compartment names
#' @param age_classes vector of upper bound for each age class
#' @param mort vector of mortaility rates for each age group
#' @param Fmat fertility matrix
#' @param vax_rates rate of vaccination
#' @param waifw who acquires infection from whom matrix 
#'              dimension: length(age_classes) x length(age_classes)
#' @param adjust_beta_flag rescale transmission rate to match homogeneous model
#'                         (using unity beta calculation)
#' @param burnin double before which negative population sizes are replenished 
#' @param print_warnings_flag print warnings if population sizes fall below 0
sir_age_structured = function(t, x, params, compartments, age_classes, mort, Fmat,
                              vax_rates, waifw, adjust_beta_flag = FALSE, 
                              print_warnings_flag = FALSE, burnin = 0){
  with(as.list(params),{
    x = x[1:(length(x)-4)] # remove beta hat variables for calculations
    if(any(x < 0)){x[which(x<0)] = 0; if(t< burnin){x[seq(2, length(x),3)] = x[seq(2, length(x),3)] + 10}; print(paste0("adj at t = ", t))} # check x is non-neg
    nage = length(age_classes)
    ncomp = length(compartments)
    aging = 1/diff(c(0,age_classes))
    x_mat = matrix(x, ncol = ncomp, nrow = nage, dimnames = list(age_classes, compartments), byrow = TRUE)
    # calculate vaccination rates (assumed to occur at birth)
    v = ifelse(is.na(vax_rates), 0, vax_rates)
    # calculate age-specific transmission rates
    beta = beta0*(beta1*cos(2*pi*(t-p))+1)
    N = sum(x_mat)
    tot_I = sum(x_mat[, "I"])
    tot_S = sum(x_mat[, "S"])
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
    der = c(dS, dI, dR)
    names(der) = sapply(compartments, function(i){paste0(i, "_", age_classes)})
    der = matrix(c(dS, dI, dR), 
                 nrow = nage, ncol = ncomp)
    der = c(t(der))
    names(der) = sapply(age_classes, function(i){paste0(compartments, "_", i)})
    der = c(der, BH = sum(lambda*x_mat[, "S"]))
    if(any(is.na(der))){browser()}
    return(list(der))
  })
}



