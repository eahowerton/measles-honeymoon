#### RUN ODE MODEL -------------------------------------------------------------
#' function to run age-specific model
#' 
#' @param age_classes vector of upper bound for each age class
#' @param compartments vector of ODE compartment names
#' @param params named vector of parameter values
#' @param mort vector of mortaility rates for each age group
#' @param fert vector of fertility rates for each age group
#' @param start_pop initial population size
#' @param vax_rates vaccination rate
#' @param waifw who acquires infection from whom matrix 
#'              dimension: length(age_classes) x length(age_classes)
#' @param IC_manual vector of proportions of the population in each compartment
#' @param max_t double max time of the simulation
#' @param dt double time step of simulation
#' @param beep_flag beep when simulation is complete
#' @param adjust_beta_flag rescale transmission rate to match homogeneous model
#'                         (using unity beta calculation)
#' @param burnin double before which negative population sizes are replenished 
#' @param print_warnings_flag print warnings if population sizes fall below 0
#' @param func function that contains age structured model
#' @return long data.frame with results for all time steps, ages, and compartments
run_ode <- function(age_classes, compartments, params, mort, fert, start_pop,
                    vax_rates = NA, waifw = NA, IC_manual = NA, max_t, dt = 1/12, 
                    beep_flag = FALSE, adjust_beta_flag = FALSE,
                    burnin = 0, print_warnings_flag = FALSE, 
                    func = sir_age_structured){
  # setup
  IC = setup_IC(start_pop = start_pop, age_classes = age_classes, 
                compartments = compartments, mort = mort, fert = fert, 
                IC_manual = IC_manual, dt = dt)
  if(any(is.na(waifw))){
    waifw = matrix(1, length(age_classes), length(age_classes)) 
  }
  # run model
  times = seq(0, max_t, dt)
  Fmat <- buildFMatrix(age.classes = age_classes, fert = fert, ncompartments = length(compartments))
  rslts <- as.data.frame(
    ode(
      y = IC,
      times = times,
      func = func,
      compartments = compartments,
      age_classes = age_classes,
      mort = mort, 
      waifw = waifw,
      vax_rates = vax_rates,
      parms = params, 
      Fmat = Fmat,
      adjust_beta_flag = adjust_beta_flag, 
      print_warnings_flag = print_warnings_flag, 
      burnin = burnin
    ))
  if(beep_flag){beep()}
  return(process_results(rslts, max_t, dt))
}

#' provide an IC vector (distributed according to stable age distribution)
#' @inheritParams run_ode
setup_IC <- function(start_pop, age_classes, compartments, mort, fert, 
                     IC_manual = NA, dt){
  # indexing - the rows for maternal, susceptible, etc
  indx_comp = rep(compartments, length(age_classes))
  IC <- rep(0,length(age_classes)*length(compartments))
  names(IC) = paste0(rep(compartments, length(age_classes)), "_", 
                     sort(rep(age_classes, length(compartments))))
  if(length(IC_manual) == length(compartments)){
    expected_stable <- findStableStruct(age.classes = age_classes, 
                                        mort = mort, fert = fert, time.step = dt)
    for(i in 1:length(compartments)){
      tmp_comp = compartments[i]
      IC[which(indx_comp == tmp_comp)] = start_pop*IC_manual[tmp_comp]*expected_stable$stable.age
    }
  }
  else if(length(IC_manual == length(IC))){
    IC = IC_manual
  }
  IC = c(IC, BH = 0)
  return(IC)
}

#' generate long data frame from ODE output
process_results <- function(rslts, max_t, dt){
  rslts_long <- rslts %>%
    mutate(BH = c(NA, diff(BH))*(1/dt)) %>%
    melt(c("time")) %>%
    tidytable::separate(variable, into = c("variable", "age"), sep = "_") %>%
    mutate(age = as.double(age))
  return(rslts_long)
}

#### HELPERS DURING SIMULATION -------------------------------------------------
#' build fertility matrix
buildFMatrix <- function(age.classes=c(1:60, seq(72,120,by=12), seq(180,600,by=60)),  
                         fert =  c(rep(0,66), rep(0.1,7)),
                         ncompartments, time.step = 1, maternal_immunity_flag = FALSE){
  nage <- length(age.classes)
  Fmat <- matrix(0,ncompartments*nage,ncompartments*nage)
  birth_compartment = 1
  for (j in 1:nage) {
    Fmat[birth_compartment,((j-1)*ncompartments+1):(j*ncompartments)] <- rep(fert[j]*time.step, ncompartments)    
  }
  return(Fmat)
}

#' find stable age structure
findStableStruct <- function(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
                             mort=c(rep(1e-9,72),1), 
                             fert =  c(rep(0,66),rep(0.1,7)), time.step = 1){
  nage <- length(age.classes)
  # aging.rate <- time.step/diff(c(0,age.classes))
  aging.rate <- 1/diff(c(0,age.classes))*time.step
  Fmat <- Tmat <- matrix(0,nage,nage)
  for (j in 1:(nage-1)) { 
    Tmat[j,j] <- (1-mort[j]*time.step)*(1-aging.rate[j])
    Tmat[j+1,j] <- (1-mort[j]*time.step)*aging.rate[j]
    Tmat[j,j] <- (1-mort[j]*time.step)*(1-aging.rate[j])
    Tmat[j+1,j] <- (1-mort[j]*time.step)*aging.rate[j]
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
    stable.age2[which(stable.age2 < 0)] = 0
    stable.age2 = stable.age2/sum(stable.age2)
    stable.age = stable.age2
  }
  lambda <- Re(eigen(Tmat+Fmat)$value[1])
  reprod.value <- Re(eigen(Tmat+Fmat)$vector[1,])
  return(list(stable.age = stable.age, lambda = lambda,
              reprod.value = reprod.value, age.classes = age.classes))
}

#### NEXT GENERATION MATRIX FUNCTIONS ------------------------------------------
#' age-specific R0/Re calculations
#' for R0 pass DFE (stable age distribution) to S
#' N is a vector of age-specific population sizes
get_Rt <- function(waifw, S, beta0, gamma, mu, N){
  S = S/N
  # if(any(length(S) != dim(waifw))){browser()}
  NGM <- beta0 / (gamma + mu) * waifw %*% diag(S)
  eigenvalues <- eigen(NGM)$values
  R0 <- max(Re(eigenvalues))
  return(R0)
}

#### FIND SCALARS --------------------------------------------------------------
#' find scalar on WAIFW matrix to achieve a given R0
find_scalar = function(s, R0, waifw, S, beta0, gamma, mu, N){
  # print(paste0("s: ", s, " R0: ", get_Rt(waifw, S, beta0*s, gamma, N)))
  diff = get_Rt(waifw, S, beta0*s, gamma, mu, N) - R0
  return(abs(diff))
}

#### WAIFWS --------------------------------------------------------------------
#' generate list of 5 sample WAIFW matrices
get_waifws = function(age_classes, background = 0.001, rescale = TRUE){
  # constant
  waifw1 <- matrix(1, length(age_classes), length(age_classes))
  # centered at age 5
  waifw2 <- get.smooth.WAIFW(age_classes = age_classes, 
                             mu = 5, sig = 0.2, gam = 0.05, delta = background)
  # cenetered at age 10
  waifw3 <- get.smooth.WAIFW(age_classes = age_classes, 
                             mu = 10, sig = 0.1, gam = 0.01, delta = background)
  # diagonal(ish)
  waifw4 <- get.smooth.WAIFW(age_classes = age_classes, 
                             mu = 7,sig = 0.6, gam = 0.05, delta = background) 
  # polymod
  waifw5 <- get.polymod.WAIFW(age_classes = age_classes, do.touch = TRUE)
  waifw = list(waifw1, waifw2, waifw3, waifw4, waifw5)
  waifw[[5]] = unname(waifw[[5]])
  if(rescale){
    # re-scale Jess's WAIFW matrices to work with my parameters
    waifw = lapply(waifw, function(i){i/mean(i)})
  }
  return(waifw)
}

#' Make a parametric smooth WAIFW
#' from Farrington, J. Amer. Stat. Assoc.; 2005, 100 p370;
#'
#' @param age_classes vector of upper bound for each age class
#' parameters defining the shape
#' @param mu age of highest contact increases with mu
#' @param gam width around equal age diagonal increases with gam
#' @param sig decreases strength in other diagonal (shrinks high trans int)
#' @param delta background homogeneous contact rate
get.smooth.WAIFW<-function(age_classes = (1:120/12),
                           mu=12.71,sig=0.69, gam=0.17, delta=0){
  n.age.cats <- length(age_classes)
  ages.to.use <- (age_classes +
                    c(0,age_classes[2:n.age.cats-1]))/2
  #gamma (p371)
  gam.func <- function(x,y,mu,sig){
    u <- (x+y)/(sqrt(2))
    vee <- 1/(sig^2)
    cval <- (sqrt(2)*mu*(1-(1/vee)))^(vee-1)
    cval <- cval*exp(1-vee)
    if (vee<1) cval <- 1
    gamma <- (1/cval)*(u^(vee-1))*exp((-vee*u)/(sqrt(2)*mu))
    return(gamma)
  }
  #b (p371); set alpha=beta here since always want symmetrical matrix
  b.func <- function(x,y,gam){
    u <- (x+y)/(sqrt(2))
    v <- (x-y)/(sqrt(2))
    alpha<-beta<-(1-gam)/(2*gam)
    b <- (((u+v)^(alpha-1))*((u-v)^(beta-1)))/(u^(alpha+beta-2))
    return(b)
  }
  #smooth value
  beta.func <- function(x,y,mu,sig, gam, delta){
    gam.func(x,y,mu,sig)*b.func(x,y,gam)+delta
  }
  betas <- outer(X=ages.to.use,Y=ages.to.use,
                 FUN=beta.func,mu=mu,sig=sig,gam=gam,delta=delta)
  return(betas)
}

#' generate POLYMOD matrix using raw contacts
#' (ignoring rescaling by population size that might be sensible)
#' 
#' @param age_classes vector of upper bound for each age class
#' @param country character of the desired country (levels possible "all", "IT"
#'                "DE" "LU" "NL" "PL" "GB" "FI" "BE"), default is "GB"
#' @param bandwidth desired smooth bandwidth
#' @param do.touch flag to include contacts involving touching
#' @param file_path path where POLYMOD data lives
get.polymod.WAIFW <- function (age_classes = (1:90),
                               country="GB",
                               bandwidth=c(3,3),
                               do.touch=FALSE, 
                               file_path = "data/polymodRaw.csv"
                               ) {
  require(KernSmooth)
  #bring in the polymod data on contacts for UK
  polymod <- read.csv(file_path)
  if (country!="all") polymod <- polymod[polymod$country==country,]
  #do touch only
  if (do.touch) polymod <- polymod[polymod[,"cnt_touch"]==1,]
  #obtain contacts, remove NAs, and make symmetrical by doubling up
  x <- cbind(polymod$participant_age, polymod$cnt_age_mean)
  x <- x[!is.na(x[,1]),]
  x <- x[!is.na(x[,2]),]
  x <- rbind(x,cbind(x[,2],x[,1]))
  #smooth monthly from 0 to 91 yrs of age and get the corresponding ages
  est<- bkde2D(x, bandwidth=bandwidth,gridsize=c(12*101, 12*101),range.x=list(c((1/24),100),c((1/24),100)))
  ages.polymod.smooth <- est$x1
  #image(ages.polymod.smooth,ages.polymod.smooth,est$fhat,xlim=c(0,10),ylim=c(0,10))
  #find which fitted ages (ages.polymod.smooth) each of the desired lower boundaries is between
  n.age.cats <- length(age_classes)
  lowpoints <-  c(0,age_classes[2:n.age.cats-1])
  index <- (findInterval(lowpoints,ages.polymod.smooth));
  #set very large ages to all be the same as the largest age
  index[index>=length(ages.polymod.smooth)]=length(ages.polymod.smooth)-1
  #extract appropriate matrix, taking smoothed estimates for the ranges
  foi.matrix <- est$fhat[index+1,index+1]
  #adjust for fact that the width of your age class should not affect the number of contacts you make
  #foi.matrix <- foi.matrix/diff(c(0,age_classes))
  colnames(foi.matrix) <- age_classes
  rownames(foi.matrix) <- age_classes
  return(foi.matrix)
}




