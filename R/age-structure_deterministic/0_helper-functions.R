#### RUN ODE MODEL -------------------------------------------------------------
#' max_t in weeks
run_ode <- function(age_classes, mort, fert, start_pop, vax_change_times = NA, vax_rates = NA, compartments,
                    waifw = NA, IC_type, IC_manual = NA, max_t, params, beep_flag = FALSE, adjust_beta_flag = FALSE, 
                    print_warnings_flag = FALSE, plot_flag = FALSE, plot_title = NA, dt = 1/12, func = sir_age_structured, burnin = 0){
  # setup
  IC = setup_IC(start_pop = start_pop, age_classes = age_classes, 
                compartments = compartments, mort = mort, fert = fert, 
                IC_type = IC_type, IC_manual = IC_manual, dt = dt)
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
      fert = fert,
      waifw = waifw,
      vax_change_times = vax_change_times, 
      vax_rates = vax_rates,
      parms = params, 
      Fmat = Fmat,
      adjust_beta_flag = adjust_beta_flag, 
      print_warnings_flag = print_warnings_flag, 
      burnin = burnin
    ))
  if(beep_flag){beep()}
  return(process_results(rslts, plot_flag, plot_title, max_t, dt))
}

# std give 2000 infected individuals across all PA and PB
setup_IC <- function(start_pop, age_classes, compartments, mort, fert, IC_type = "std", IC_manual = NA, dt, rel_size = NA){
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
    IC[which(indx_comp == "S1")] = start_pop*rel_size[1]*0.059/length(age_classes) 			 # HACK FOR NOW... FIX THIS...
    IC[which(indx_comp == "I1")] = 0.001*rel_size[1]/length(age_classes)
    IC[which(indx_comp == "R1")] = start_pop*rel_size[1]*0.94/length(age_classes)
    IC[which(indx_comp == "S2")] = start_pop*rel_size[2]*0.059/length(age_classes) 			 # HACK FOR NOW... FIX THIS...
    IC[which(indx_comp == "I2")] = 0.001*rel_size[2]/length(age_classes)
    IC[which(indx_comp == "R2")] = start_pop*0.94*rel_size[2]/length(age_classes)
  }
  if(IC_type == "std-noI"){
    # setup initial conditions
    IC[which(indx_comp == "S")] = start_pop*0.06/length(age_classes) 			 # susceptibles
    IC[which(indx_comp == "R")] = start_pop*0.94/length(age_classes)
    IC[which(indx_comp == "S1")] = start_pop*rel_size[1]*0.06/length(age_classes) 			 # HACK FOR NOW... FIX THIS...
    IC[which(indx_comp == "R1")] = start_pop*rel_size[1]*0.94/length(age_classes)
    IC[which(indx_comp == "S2")] = start_pop*rel_size[2]*0.06/length(age_classes) 			 # HACK FOR NOW... FIX THIS...
    IC[which(indx_comp == "R2")] = start_pop*rel_size[2]*0.94/length(age_classes)
  }
  else if(IC_type == "stable-age"){
    expected_stable <- findStableStruct(age.classes = age_classes, mort = mort, fert = fert, time.step = dt)
    IC[which(indx_comp == "S")] = start_pop*0.059*expected_stable$stable.age 			 # susceptibles
    IC[which(indx_comp == "I")] = start_pop*0.001*expected_stable$stable.age 			 # infecteds
    IC[which(indx_comp == "R")] = start_pop*0.94*expected_stable$stable.age 			 # recovereds
    IC[which(indx_comp == "S1")] = start_pop*rel_size[1]*0.059*expected_stable$stable.age 			 # susceptibles
    IC[which(indx_comp == "I1")] = start_pop*rel_size[1]*0.001*expected_stable$stable.age 			 # infecteds
    IC[which(indx_comp == "R1")] = start_pop*rel_size[1]*0.94*expected_stable$stable.age 			   # recovereds
    IC[which(indx_comp == "S2")] = start_pop*rel_size[2]*0.059*expected_stable$stable.age 			 # susceptibles
    IC[which(indx_comp == "I2")] = start_pop*rel_size[2]*0.001*expected_stable$stable.age 			 # infecteds
    IC[which(indx_comp == "R2")] = start_pop*rel_size[2]*0.94*expected_stable$stable.age 			   # recovereds
    if(any(IC < 0)){browser()}
    if(abs(sum(IC) - start_pop) > 1e-6){browser()}
  }
  else if(IC_type == "stable-age-noI"){
    expected_stable <- findStableStruct(age.classes = age_classes, mort = mort, fert = fert, time.step = dt)
    IC[which(indx_comp == "S")] = start_pop*0.06*expected_stable$stable.age 			 # susceptibles
    IC[which(indx_comp == "R")] = start_pop*0.94*expected_stable$stable.age 			 # recovereds
    IC[which(indx_comp == "S1")] = start_pop*0.06*rel_size[1]*expected_stable$stable.age 			 # susceptibles
    IC[which(indx_comp == "R1")] = start_pop*0.94*rel_size[1]*expected_stable$stable.age 			   # recovereds
    IC[which(indx_comp == "S2")] = start_pop*0.06*rel_size[2]*expected_stable$stable.age 			 # susceptibles
    IC[which(indx_comp == "R2")] = start_pop*0.94*rel_size[2]*expected_stable$stable.age 			   # recovereds
    if(any(IC < 0)){browser()}
    if(abs(sum(IC) - start_pop) > 1e-6){browser()}
  }
  # otherwise provide an IC vector (distributed according to stable age distribution)
  else if(IC_type == "manual")
  {
    if(length(IC_manual) == length(compartments)){
      expected_stable <- findStableStruct(age.classes = age_classes, mort = mort, fert = fert, time.step = dt)
      for(i in 1:length(compartments)){
        tmp_comp = compartments[i]
        IC[which(indx_comp == tmp_comp)] = start_pop*IC_manual[tmp_comp]*expected_stable$stable.age
      }
    }
    else if(length(IC_manual == length(IC))){
      IC = IC_manual
    }
  }
  IC = c(IC, C = 0, BH = 0, BHs = 0, BHi = 0, BHb = 0)
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

#### HELPERS DURING SIMULATION -------------------------------------------------
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
    stable.age2[which(stable.age2 < 0)] = 0
    stable.age2 = stable.age2/sum(stable.age2)
    stable.age = stable.age2
  }
  lambda <- Re(eigen(Tmat+Fmat)$value[1])
  reprod.value <- Re(eigen(Tmat+Fmat)$vector[1,])
  return(list(stable.age = stable.age, lambda = lambda,
              reprod.value = reprod.value, age.classes = age.classes))
}

# for R0 pass DFE (stable age distribution) to S
# N is a vector of age-specific population sizes
get_Rt <- function(waifw, S, beta0, gamma, N){
  S = S/N
  NGM <- beta0 / gamma * waifw %*% diag(S)
  eigenvalues <- eigen(NGM)$values
  R0 <- max(Re(eigenvalues))
  return(R0)
}

compute_extinction_prob <- function(waifw, S, beta0, gamma, N, age_classes,
                                    tol = 1e-8, max_iter = 1000) {
  S = S/N
  NGM <- beta0 / gamma * waifw %*% diag(S)
  q <- rep(0.5, nrow(NGM))  # start at full extinction
  for (i in 1:max_iter) {
    q_new <- exp(NGM %*% (q - 1))
    if (max(abs(q_new - q)) < tol) break
    q <- q_new
  }
  return(data.frame(age = age_classes, extinction_prob = q, outbreak_prob = 1 - q))
}

#### FIND SCALARS --------------------------------------------------------------
match_on_age_or_inc = function(age_classes, mort, fert, start_pop, compartments, 
                               vax_pct = 0, params, waifws, max_t = 100, 
                               optim_quantity = "tot_I"){
  IC = setup_IC(start_pop, age_classes, compartments, mort, fert, IC_type = "stable-age")
  Fmat <- buildFMatrix(age.classes = age_classes, fert = fert, ncompartments = length(compartments))
  # run ODE with flat WAIFW
  flat = runsteady(
    y = IC, times = c(0, 500), func = sir_age_structured, parms = params, # runsteady arguments
    compartments = compartments, age_classes = age_classes, mort = mort, fert = fert,
    vax_change_times = c(0), vax_rates = c(vax_pct),
    waifw = matrix(1, length(age_classes), length(age_classes)),
    Fmat = Fmat, adjust_beta_flag = FALSE, print_warnings_flag = FALSE
  )
  which_I = which(compartments == "I")
  # get age distribution of cases
  age_dist_target = flat$y[seq(which_I, length(flat$y)-4, length(compartments))] # assuming I = 3rd compar
  # loop through beta multipliers for all other waifw 
  waifw_scalars = vector("list", length(waifws))
  for(i in 1:length(waifws)){
    print(paste0("starting waifw ", i, "/", length(waifws)))
    tmp =  tryCatch(
      optimize(f = match_one_waifw, tol = 1e-1, interval = c(0, 3), # optim arguments  lower = 0, 
               IC = IC, Fmat = Fmat, age_classes = age_classes, mort = mort, fert = fert, 
               compartments = compartments, params = params,
               vax_pct = vax_pct, max_t = max_t,
               age_dist_flat = age_dist_target, new_waifw = waifws[[i]],
               optim_quantity = optim_quantity)
      , error = function(msg){
        return(NA)#data.frame(minimum = NA, diff = NA))
      })
    if(any(is.na(tmp))){
      waifw_scalars[[i]] = data.frame(best_scalar = NA, 
                                      diff = NA)
    }
    else{
      waifw_scalars[[i]] = data.frame(best_scalar = tmp$minimum, 
                                      diff = tmp$objective)
    }
  }
  return(bind_rows(waifw_scalars, .id = "waifw_id") %>% mutate(waifw_id = as.integer(waifw_id) + 1))
}

match_one_waifw = function(beta_scalar, age_classes, mort, fert, IC, Fmat,
                           compartments, params, age_dist_flat, new_waifw, 
                           vax_pct, max_t, optim_quantity = "tot_I"){
  tmp_pars = params
  tmp_pars["beta0"] = tmp_pars["beta0"]*beta_scalar
  tmp_out = tryCatch(
    #this is the chunk of code we want to run
    {runsteady(y = IC, times = c(0, 500), func = sir_age_structured, parms = tmp_pars, # runsteady arguments
               compartments = compartments, age_classes = age_classes, mort = mort, fert = fert,
               waifw = new_waifw, vax_change_times = c(0), vax_rates = c(vax_pct),
               Fmat = Fmat, adjust_beta_flag = FALSE, print_warnings_flag = FALSE
    )
      #when it throws an error, the following block catches the error
    }, error = function(msg){
      return(data.frame(time = NA))
    })
  if(length(unlist(tmp_out)) == 1){
    print("Numerical error, trying again")
    tmp_pars = params
    tmp_pars["beta0"] = tmp_pars["beta0"]*(beta_scalar + 1e-4)
    tmp_out = tryCatch(
      #this is the chunk of code we want to run
      {runsteady(y = IC, times = c(0, 500), func = sir_age_structured, parms = tmp_pars, # runsteady arguments
                 compartments = compartments, age_classes = age_classes, mort = mort, fert = fert,
                 waifw = new_waifw, vax_change_times = c(0), vax_rates = c(vax_pct),
                 Fmat = Fmat, adjust_beta_flag = FALSE, print_warnings_flag = FALSE
      )
        #when it throws an error, the following block catches the error
      }, error = function(msg){
        return(data.frame(minimum = NA, diff = NA))
      })
    if(length(unlist(tmp_out))== 1){return(NA)}
  }
  which_I = which(compartments == "I")
  age_dist_new_waifw = tmp_out$y[seq(which_I, length(tmp_out$y)-4, length(compartments))] # assuming I = 3rd compar
  if(any(age_dist_new_waifw < -1e-2)){print("negative values, trying again");return(1e6*beta_scalar)}
  if(optim_quantity == "tot_I"){
    print(paste0("scalar: ", round(beta_scalar,3), " diff: ", round(abs(sum(age_dist_new_waifw) - sum(age_dist_flat)),3)))
    return(abs(sum(age_dist_new_waifw) - sum(age_dist_flat))) # minimize absolute value of difference
  }
  else if(optim_quantity == "mean_age"){
    bin_width = diff(c(0, age_classes))
    flat_mean_age = sum(age_dist_flat*(age_classes-bin_width/2))/sum(age_dist_flat)
    return(abs(new_waifw_mean_age - flat_mean_age))
  }
}

#### WAIFWS --------------------------------------------------------------------
get_waifws = function(age_classes, background = 0.001, rescale = TRUE){
  # constant
  waifw1 <- matrix(1, length(age_classes), length(age_classes))
  # centered at age 5
  waifw2 <- get.smooth.WAIFW(age.class.boundaries = age_classes, 
                             mu = 5, sig = 0.2, gam = 0.05, delta = background)
  # cenetered at age 10
  waifw3 <- get.smooth.WAIFW(age.class.boundaries = age_classes, 
                             mu = 10, sig = 0.1, gam = 0.01, delta = background)
  # diagonal(ish)
  waifw4 <- get.smooth.WAIFW(age.class.boundaries = age_classes, 
                             mu = 7,sig = 0.6, gam = 0.05, delta = background) 
  # polymod
  waifw5 <- get.polymod.WAIFW(age.class.boundaries = age_classes, do.touch = TRUE)
  waifw = list(waifw1, waifw2, waifw3, waifw4, waifw5)
  waifw[[5]] = unname(waifw[[5]])
  if(rescale){
    # re-scale Jess's WAIFW matrices to work with my parameters
    waifw = lapply(waifw, function(i){i/mean(i)})
  }
  return(waifw)
}

# Make a parametric smooth WAIFW
# from Farrington, J. Amer. Stat. Assoc.; 2005, 100 p370;
#
#  parameters -
#     age class boundries - the upper age limit for each age class in years,
#     parameters defining the shape:
#        mu: age of highest contact increases with mu
#        gam: width around equal age diagonal increases with gam
#        sig: decreases strength in other diagonal (shrinks high trans int)
#        delta:  background homogeneous contact rate
#
#Returns -
#   a smooth WAIFW matrix
get.smooth.WAIFW<-function(age.class.boundaries = (1:120/12),
                           mu=12.71,sig=0.69, gam=0.17, delta=0){
  n.age.cats <- length(age.class.boundaries)
  ages.to.use <- (age.class.boundaries +
                    c(0,age.class.boundaries[2:n.age.cats-1]))/2
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

#Make a WAIFW matrix based on Polymod (ignoring rescaling by population
#size that might be sensible - just use raw contacts)
#
#Parameters -
#   age class boundries - the upper age limit for each age class in YEARS
#   the desired country - levels possible "all", "IT" "DE" "LU" "NL" "PL" "GB" "FI" "BE",
#         default is "GB"
#   bandwidth - desired smooth bandwidth - default=c(3,3)
#   do.touch - do you want to only include contacts involving touching? default is FALSE
#Returns -
#   a WAIFW matrix based on the Polymod results from chosen location with row and col
#   names indicating age classes

get.polymod.WAIFW <- function (age.class.boundaries = (1:90),
                               country="GB",
                               bandwidth=c(3,3),
                               do.touch=FALSE) {
  require(KernSmooth)
  
  #bring in the polymod data on contacts for UK
  polymod <- read.csv("data/polymodRaw.csv")
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
  n.age.cats <- length(age.class.boundaries)
  lowpoints <-  c(0,age.class.boundaries[2:n.age.cats-1])
  index <- (findInterval(lowpoints,ages.polymod.smooth));
  #set very large ages to all be the same as the largest age
  index[index>=length(ages.polymod.smooth)]=length(ages.polymod.smooth)-1
  
  #extract appropriate matrix, taking smoothed estimates for the ranges
  foi.matrix <- est$fhat[index+1,index+1]
  
  #adjust for fact that the width of your age class should not affect the number of contacts you make
  #foi.matrix <- foi.matrix/diff(c(0,age.class.boundaries))
  
  colnames(foi.matrix) <- age.class.boundaries
  rownames(foi.matrix) <- age.class.boundaries
  return(foi.matrix)
}




