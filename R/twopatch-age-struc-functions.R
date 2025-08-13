sirtwopatch_age_structured <- function(t, x, parms, compartments, age_classes, mort, vax_change_times, vax_rates, Fmat,
                               fert, waifw, adjust_beta_flag = FALSE, print_warnings_flag = FALSE){
  with(as.list(parms),{
    x = x[1:(length(x)-5)] # remove beta hat variables for calculations
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
    N1 = sum(x_mat[, "S1"] + x_mat[, "I1"] + x_mat[, "R1"])
    N2 = sum(x_mat[, "S2"] + x_mat[, "I2"] + x_mat[, "R2"])
    # N1 = N1 + pi*N2
    # N2 = pi*N1 + N2
    N1 = ifelse(N1 == 0, 1, N1)
    N2 = ifelse(N2 == 0, 1, N2)
    lambda1 <- waifw%*%(beta/N1*(x_mat[, "I1"] + pi*x_mat[, "I2"]))
    lambda2 <- waifw%*%(beta/N2*(pi*x_mat[, "I1"] + x_mat[, "I2"]))
    # fertility
    N_fert <- Fmat %*% x # assume all births go into S1
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
    dS1 <- - lambda1 * x_mat[, "S1"] - delta * x_mat[, "S1"]  - mort * x_mat[, "S1"] +
      N_fert[, "S1"] * v * (1 - ve) + N_age_in[, "S1"] - N_age_out[, "S1"]
    dI1 <- lambda1 * x_mat[, "S1"] - gamma * x_mat[, "I1"] + delta * x_mat[, "S1"] - mort * x_mat[, "I1"] + 
      N_fert[, "I1"] + N_age_in[, "I1"] - N_age_out[, "I1"]
    dR1 <- gamma * x_mat[, "I1"] - mort * x_mat[, "R1"] +
      N_fert[, "R1"] + N_fert[, "S1"] * v * ve + N_age_in[, "R1"] - N_age_out[, "R1"]
    dS2 <- - lambda2 * x_mat[, "S2"] - delta * x_mat[, "S2"] - mort * x_mat[, "S2"] +
      N_fert[, "S1"] * (1 - v) + N_age_in[, "S2"] - N_age_out[, "S2"]
    dI2 <- lambda2 * x_mat[, "S2"] - gamma * x_mat[, "I2"] + delta * x_mat[, "S2"] - mort * x_mat[, "I2"] +
      N_fert[, "I2"] + N_age_in[, "I2"] - N_age_out[, "I2"]
    dR2 <- gamma * x_mat[, "I2"] - mort * x_mat[, "R2"] +
      N_fert[, "R2"] + N_age_in[, "R2"] - N_age_out[, "R2"] 
    der <- c(dS1, dI1, dR1, dS2, dI2, dR2)
    names(der) <- sapply(compartments, function(i){paste0(i, "_", age_classes)})
    der <- matrix(c(dS1, dI1, dR1, dS2, dI2, dR2), 
                  nrow = nage, ncol = ncomp)
    der <- c(t(der))
    names(der) <- sapply(age_classes, function(i){paste0(compartments, "_", i)})
    der <- c(der, C = sum(lambda1*x_mat[, "S1"] + lambda2*x_mat[, "S2"]), 
             BH = 0, BHi = 0, BHs = 0, BHb = 0)
    return(list(der))
  })
}

# N1/2 are vectors of age and patch specific population sizes
get_Rt_twopatch <- function(waifw, S1, S2, N1, N2, beta0, gamma, pi){
  if(any(c(N1, N2) == 0)){return(1e6)}
  S = c(S1/sum(N1), S2/sum(N2))
  waifw_expand = rbind(
    cbind(waifw, pi*waifw), # top row
    cbind(pi*waifw, waifw)  # bottom row
  )
  NGM <- beta0 / gamma * waifw_expand %*% diag(S)
  eigenvalues <- eigen(NGM)$values
  R0 <- max(Re(eigenvalues))
  return(R0)
}


get_Rt_twopatch <- function(waifw, S1, S2, N1, N2, beta0, gamma, pi){
  if(any(c(N1, N2) == 0)){return(1e6)}
  S = c(S1/(sum(N1) + pi*sum(N2)), S2/(pi*sum(N1) + sum(N2)))
  waifw_expand = rbind(
    cbind(waifw, pi*waifw), # top row
    cbind(pi*waifw, waifw)  # bottom row
  )
  NGM <- beta0 / gamma * waifw_expand %*% diag(S)
  eigenvalues <- eigen(NGM)$values
  R0 <- max(Re(eigenvalues))
  return(R0)
}

