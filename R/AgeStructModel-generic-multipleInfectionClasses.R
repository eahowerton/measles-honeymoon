
## In this code, we distinguish between primary, secondary and tertiary or higher infections 
## This more closely maps to EVD68, and to RSV ##

#setwd("/Users/cmetcalf/Dropbox (Princeton)/RSV-Age-structure/")

## Function to bring in this chunk of code. 
#
reload.source <- function(fname=""){
    source(paste(fname,"source/AgeStructModel-generic-multipleInfectionClasses.R",sep=""))
}



## Function to build survival, aging and infectious transitions (T matrix called 'U' in the methods)
## In this instance, the epidemiological categories are: Maternaly Immune, Suscepticible, Infected via a live vaccine, Susceptbile following live, 
##								Infected_1, Susceptible_1,Infected_2, Susceptible_2, Infected_3, Susceptible_3, Infected_4, Vaccinated 
#
buildTMatrix <- function(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)),  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years (could extend)
                                    mort=c(rep(1e-9,72),1),						#neglgible death till age 50 (could make more reaslitic)
                                    beta=c(0.1),							#transmission (can be single value or vector of same length as age)
                                    waning.maternal.par=0.5,					#rate of waning of maternal immunity
                                    recovery.infection=c(1,1,1,1), 					#rate of recovery frmo infection for each of the four classes
				  vaccination=c(rep(0,68),rep(0.5,5)),				#vector of length ages (e.g., vaccinate at age 20, vaccinate 50% of the population)		
				  live.vaccination=rep(0,73),					#vector of length ages to capture vaccination of young kids witha live vaccine		
                                    waning.live.vacc=0.2, 					        #rate of waning of life vaccine
                                    waning.other.vacc=0.2, 					        #rate of waning of other vaccine
				  effects.on.trans=c(1,1,1,1),					#how much transmission attenuated for subsequent infectoins
				  time.step=1/4){						#time-step is 1/4 of a month, aka 1 week (note slight weirdness with month/year 12 months = 52 weeks)

    ## matrix dimensions, aging rate 
    nage <- length(age.classes)
    aging.rate <- time.step/diff(c(0,age.classes))

    ## adjustments for cases where beta is constant across age and only one number provided
    if (length(beta)==1) beta <- rep(beta,nage)

    ## adjustment for cases where only one waning of immunity is provided - there should be one for each type of ifnection
    if (length(recovery.infection)==1) recovery.infection <- rep(recovery.infection,4)
    if (length(effects.on.trans)==1) effects.on.trans <- rep(effects.on.trans,4)

    ## waning of maternal immunity: first calculate the overall pattern, i.e., proportion of individuals having 'survived' waning in every age class	
    pmaternal <- exp(-c(0,age.classes)*waning.maternal.par)

    ## What leaks out specifically in each age class - taking bottom and top of age class and getting proportion that have gone
    waning.maternal <- (pmaternal[1:(length(pmaternal)-1)]-pmaternal[2:length(pmaternal)])/pmaternal[1:(length(pmaternal)-1)]	

    #two matrices for storage	    
    mat1 <- matrix(0,12,12)
    Tmat <- matrix(0,12*nage,12*nage)

   

    #loop over ages
    for (j in 1:nage) {
        #fill in epi matrix	
        mat1[] <- 0

        #moving out of maternal immunity 
        mat1[1,1] <- (1-waning.maternal[j])*(1-live.vaccination[j])
        mat1[2,1] <- waning.maternal[j]*(1-live.vaccination[j])
        mat1[3,1] <- live.vaccination[j]*(1-waning.maternal[j])

	#moving out of susceptible
        mat1[2,2] <- (1-beta[j])*(1-vaccination[j])*(1-live.vaccination[j]) #stay
        mat1[3,2] <- (1-beta[j])*(1-vaccination[j])*live.vaccination[j]	#get the live vaccine
        mat1[5,2] <- beta[j]*(1-vaccination[j])*(1-live.vaccination[j])	#get infected
        mat1[12,2] <- vaccination[j]*(1-live.vaccination[j])	#get the later life vaccine (which shapes maternal imm)

	#moving out of vaccinated with live vaccine
        mat1[3,3] <- (1-waning.live.vacc)
        mat1[4,3] <- waning.live.vacc

	#Moving out of susceptible following vaccination with live vaccine
        mat1[4,4] <- (1-beta[j]*effects.on.trans[1])*(1-vaccination[j])
        mat1[7,4] <- beta[j]*effects.on.trans[1]*(1-vaccination[j])
        mat1[12,4] <- vaccination[j]

	#Moving out of infected_1 (assume don't vaccinate infected)
	mat1[5,5] <- (1-recovery.infection[1])
	mat1[6,5] <- recovery.infection[1]

	#Moving out of susceptible_1
        mat1[6,6] <- (1-beta[j]*effects.on.trans[2])*(1-vaccination[j]) #stay
        mat1[7,6] <- beta[j]*effects.on.trans[2]*(1-vaccination[j])	#get infected
        mat1[12,6] <- vaccination[j]	#get the later life vaccine (which shapes maternal imm)

      #  print(j)
	#print(sum(mat1[,6]))

	#Moving out of infected_2
	mat1[7,7] <- (1-recovery.infection[2])
	mat1[8,7] <- recovery.infection[2]

	#Moving out of susceptible_2
        mat1[8,8] <- (1-beta[j]*effects.on.trans[3])*(1-vaccination[j]) #stay
        mat1[9,8] <- beta[j]*effects.on.trans[3]*(1-vaccination[j])	#get infected
        mat1[12,8] <- vaccination[j]	#get the later life vaccine (which shapes maternal imm)

	#Moving out of infected_3
	mat1[9,9] <- (1-recovery.infection[3])
	mat1[10,9] <- recovery.infection[3]

	#Moving out of susceptible_3
        mat1[10,10] <- (1-beta[j]*effects.on.trans[4])*(1-vaccination[j]) #stay
        mat1[11,10] <- beta[j]*effects.on.trans[4]*(1-vaccination[j])	#get infected
        mat1[12,10] <- vaccination[j]	#get the later life vaccine (which shapes maternal imm)

	#Moving out of infected_4
	mat1[11,11] <- (1-recovery.infection[4])
	mat1[10,11] <- recovery.infection[4]

	#Moving out of vaccinated
	mat1[12,12] <- (1-waning.other.vacc)
	mat1[10,12] <- waning.other.vacc


        #put in surv
        surv <- rep(1-mort[j],12);
   
        #fill in Tmatrix
        if (j!=nage) {
            #go into next age class
            Tmat[(j*12+1):(j*12+12),((j-1)*12+1):(j*12)] <- mat1*surv*aging.rate[j]

            #stay in this age class - only lose maternal immunity when change an age class
            mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
            Tmat[((j-1)*12+1):(j*12),((j-1)*12+1):(j*12)] <- mat2*surv*(1-aging.rate[j])

        } else {
            #stay in the last age class
            Tmat[((j-1)*12+1):(j*12),((j-1)*12+1):(j*12)] <- mat1*surv
        }}

	return(Tmat)

}


## Funtion to build fertility transitions
#
buildFMatrix <- function(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)),  
                         fert =  c(rep(0,66),rep(0.1,7))){ 	#one for every age class. start reproducing at age 20, and here assume ~ flat
	nage <- length(age.classes)

	Fmat <- matrix(0,12*nage,12*nage)
	for (j in 1:nage) {
		#Fmat[1,((j-1)*12+1):(j*12)] <- c(0,0,fert[j],0,fert[j],0,fert[j],0,fert[j],0,fert[j],fert[j]) 	#the maternally immune born from infected andvacc
		#Fmat[2,((j-1)*12+1):(j*12)]<- c(fert[j],fert[j],0,fert[j],0,fert[j],0,fert[j],0,fert[j],0,0)  	#the susceptible (mom didn't get sick, kid born susceptible)

		Fmat[1,((j-1)*12+1):(j*12)] <- c(rep(0,11),fert[j]) 					#option - only vaccinated have maternal immunity 
		Fmat[2,((j-1)*12+1):(j*12)] <- c(rep(fert[j],11),0) 					#everyone else has susceptible babies

	}


	return(Fmat)
}


## Function to find stable age structure; population rate of increase, etc; from fertility and mortality profiles (no classification by epidemiological class)

findStableStruct <- function(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
			mort=c(rep(1e-9,72),1), 
                         fert =  c(rep(0,66),rep(0.1,7)), time.step=1/4){


	   nage <- length(age.classes)
  	   aging.rate <- time.step/diff(c(0,age.classes))

	   Fmat <- Tmat <- matrix(0,nage,nage)
	   for (j in 1:(nage-1)) { 
		Tmat[j,j] <- (1-mort[j])*(1-aging.rate[j])
		Tmat[j+1,j] <- (1-mort[j])*aging.rate[j]
	   }
	   j <- nage	
	   Tmat[j,j] <- (1-mort[j])
	   Fmat[1,] <- fert	

	   stable.age <- Re(eigen(Tmat+Fmat)$vector[,1])
	   stable.age <- stable.age/sum(stable.age)
	   lambda <- Re(eigen(Tmat+Fmat)$value[1])
	   reprod.value <- Re(eigen(Tmat+Fmat)$vector[1,])

	return(list(stable.age=stable.age,lambda=lambda,reprod.value=reprod.value,age.classes=age.classes))
}

## Function to adjust fertility until you get a desired number of babies produced at a focal popsize at the equilbrium age structure
##
findFert <- function(target.births,
			fert.multiplier=seq(1,100,length=10),  ## adjust this range based on initial estimates, make finer scale fo rmore precision, etc
			popsize,
			age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
			mort=c(rep(1e-9,72),1), 
                         fert =  c(rep(0,66),rep(0.1,7)), 
			time.step=1/4){

	    babies <- rep(NA,length(fert.multiplier))

	    for (j in 1:length(fert.multiplier)) { 

	    stable.struct <- findStableStruct(age.classes=age.classes, 
				mort=mort,fert = fert*fert.multiplier[j], time.step=time.step)
	    babies[j] <- sum(stable.struct$stable.age*popsize*fert*fert.multiplier[j])

	    }

	plot(fert.multiplier,babies, xlab="Multiplier", ylab="Total babies")
	abline(h=target.births)
	best <- which(abs(babies-target.births)==min(abs(babies-target.births)))
	abline(v=fert.multiplier[best])
	
	stable.struct <- findStableStruct(age.classes=age.classes, 
				mort=mort,fert = fert*fert.multiplier[best], time.step=time.step)

	print(stable.struct$lambda) ##just out of curiosity - is the popuation growing or shrinking? 

	return(fert.multiplier[best])

}




# Function to iterate out to equilibrium
#
#
iterateSIR <- function(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)),  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                        mort=c(rep(0,72),1),					#neglgible death till age 50 (could make more reaslitic)
                        waifw=matrix(0.0002,73,73),					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
			alpha=0.2,							#parameter governing sine wave seasonality in transmission	
			discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                         waning.maternal.par=0.5,						#rate of waning of maternal immunity
                         recovery.infection=c(1,1,1,1), 					#rate of recovery frmo infection for each of the four classes
			vaccination=c(rep(0,68),rep(0,5)),				#vector of length ages (e.g., vaccinate at age 20, vaccinate 50% of the population)		
			live.vaccination=rep(0,73),					#vector of length ages to capture vaccination of young kids witha live vaccine		
                          waning.live.vacc=0.2, 					        #rate of waning of life vaccine
                          waning.other.vacc=0.2, 					        #rate of waning of other vaccine
			effects.on.trans=c(1,1,1,1),					#how much transmission attenuated for susceptibles who have previously been infected
			contribution.to.trans=c(1,0.5,0.1,0.1),					#how much transmission is reduced from primary, secondary, tertiary, etc infections
			time.step=1/4,							#time-step (1/4 of a month = 1 week)
                         fert= c(rep(0,66),rep(0.0008,7)),					#pattern of fertility over age
			alpha.births = 0,						#pattern of seasonality in births	
			intervention=c(250,250+52*2,1),					#start and end time point and reduction in transmission of the NPIs
			burnin=10,Tmax=150, 						#control variables on time
			do.plot=FALSE,							#plot or not
			start.pop=c(),							#starting population structure and #could impose, otherwise takes the eq structure 
			start.pop.size=100000,						#starting total population size
			gamma=0.97){							#this parameter mitigates explosiveness emerging from discrete time. 


    #indexing - the rows for maternal, susceptible, etc
    m.inds <- seq(1,length(age.classes)*12,by=12)
    s.inds <- m.inds+1
    i.indsv <- m.inds+2 	#live vaccined
    s.indsv <- m.inds+3	#recover from live vaccine
    i.inds1 <- m.inds+4	#first infection
    s.inds1 <- m.inds+5	#recover from 1st
    i.inds2 <- m.inds+6  #second infection...
    s.inds2 <- m.inds+7
    i.inds3 <- m.inds+8
    s.inds3 <- m.inds+9
    i.inds4 <- m.inds+10
    v.inds <- m.inds+11  #other vaccine


    #get Fmatrix (basic won't change - might multiply by seasonal flucutation, but can happen as single parameter)
    Fmat <- buildFMatrix(age.classes=age.classes,fert=fert)

    #initiate populations
    find.expected.stable <- findStableStruct(age.classes=age.classes,mort= mort,fert =  fert, time.step=time.step)
    if (length(start.pop)==0) start.pop <- start.pop.size*find.expected.stable$stable.age

    nt<-rep(0,length(age.classes)*12)
    nt[s.inds] <- start.pop  			 #susceptibles
    nt[i.inds1[1:60]] <- 10 			 #seed some infected (otherwise infection can't take off)
    nt[s.inds[1:60]] <- nt[s.inds[1:60]]-10 	 #remove to keep pop size the same
    Rt<-store.beta <- rep(NA,Tmax)		 #storage
    nt.store <- matrix(NA,length(nt),Tmax)

    #seasaonl transmission - each has a mean / median of 1
    if (length(discrete.seas)==0) seas <- (1+alpha*cos(2*pi*(1:52)/52))	else { seas <- discrete.seas/median(discrete.seas); beta <-1 }

    # modulation of magnitude of transmission (i.e., NPIs)
    transmission.reduction <- rep(1,Tmax)
    transmission.reduction[intervention[1]:intervention[2]] <- intervention[3]
 

    for (t in 1:Tmax) {

	#extract rate of transmission to everage age
	phi <- (1-exp(-transmission.reduction[t]*seas[1+t%%52]*waifw%*%((contribution.to.trans[1]*nt[i.inds1]+contribution.to.trans[2]*nt[i.inds2]+
					contribution.to.trans[3]*nt[i.inds3]+contribution.to.trans[4]*nt[i.inds4])^gamma))); 
	beta.trans <- t(phi); 

	#print(range(beta.trans))

        #build the T matrix 
        Tmat <- buildTMatrix(age.classes=age.classes,
                                        mort=mort,   	
                                        beta=beta.trans,
                                        waning.maternal.par=waning.maternal.par,
                                        recovery.infection=recovery.infection, 	
				      vaccination=vaccination,				#vector of length ages (e.g., vaccinate at age 20, vaccinate 50% of the population)		
				      live.vaccination=live.vaccination,					#vector of length ages to capture vaccination of young kids witha live vaccine		
                                        waning.live.vacc=waning.live.vacc, 					        #rate of waning of life vaccine
                                        waning.other.vacc=waning.other.vacc, 					        #rate of waning of other vaccine
				      effects.on.trans=effects.on.trans,					#how much transmission attenuated for subsequent infectoins
                                        time.step=time.step)

	

	#introduce seasonality in births (if alpha.births=0, this will just be x1)
	seas.births <- (1+alpha.births*cos(2*pi*t/52))	

	#use matrix multiplication to get the population structure at the next time-step
	nt1<-(Tmat+Fmat*seas.births) %*% nt
	
	#print(nt1[1:14])

        sum.nt1<-sum(nt1)
        Rt[t]<-log(sum.nt1)	#pop growth rate
        nt<-nt1 			#update the population 
        nt.store[,t] <- nt	#store pop
        store.beta[t] <- max(beta.trans)	#stoer transmissino

        #if earlier than burnin, and no infecteds, introduce some
        if (sum(nt[i.inds1])<1 & t<burnin) {nt[i.inds1] <- 1/length(age.classes);nt[s.inds] <- pmax(nt[s.inds]-1/length(age.classes),0) }
    }
 
    if (do.plot) { ## mostly sanity checks / diagnoses
        par(mfrow=c(2,1))
        image((burnin:Tmax)/52,age.classes/12,log(t(nt.store[s.inds,burnin:Tmax])+1), xlab="Time", ylab="Age", main="susceptible", ylim=c(0,15))
        contour((burnin:Tmax)/52,age.classes/12,t(nt.store[s.inds,burnin:Tmax])+1,add=TRUE)
        image((burnin:Tmax)/52,age.classes/12,log(t(nt.store[i.inds1,burnin:Tmax])+1), xlab="Time", ylab="Age", main="infected", ylim=c(0,15))
        contour((burnin:Tmax)/52,age.classes/12,t(nt.store[i.inds1,burnin:Tmax])+1,add=TRUE)


	par(mfrow=c(2,2))
  	plot((burnin:Tmax)/52, colSums(nt.store[i.inds1,burnin:Tmax]), type="l",col=1,xlab="Time", ylab="Total cases")
	abline(v=c(0:100),lty=3)
	plot((burnin:Tmax)/52, colSums(nt.store[m.inds,burnin:Tmax])+colSums(nt.store[m.inds+1,burnin:Tmax])+
				colSums(nt.store[m.inds+2,burnin:Tmax])+colSums(nt.store[m.inds+3,burnin:Tmax])+
				colSums(nt.store[m.inds+4,burnin:Tmax])+colSums(nt.store[m.inds+5,burnin:Tmax])+
				colSums(nt.store[m.inds+6,burnin:Tmax])+colSums(nt.store[m.inds+7,burnin:Tmax])+
				colSums(nt.store[m.inds+8,burnin:Tmax])+colSums(nt.store[m.inds+9,burnin:Tmax])+
				colSums(nt.store[m.inds+10,burnin:Tmax])+colSums(nt.store[m.inds+11,burnin:Tmax]), 
				type="l",col=1,xlab="Time", ylab="Total pop size")#, xlim=c(120,134))
	
	plot(age.classes/12,nt.store[m.inds,Tmax]/diff(c(0,age.classes)), type="l", xlab="Ages", ylab="Number", xlim=c(0,10))
	points(age.classes/12,nt.store[m.inds,burnin]/diff(c(0,age.classes)), type="l")

	points(age.classes/12,nt.store[i.inds1,Tmax]/diff(c(0,age.classes)),type="l",col=2)
	points(age.classes/12,nt.store[i.inds1,burnin]/diff(c(0,age.classes)),col=2, type="l")
	legend("topright",legend=c("Infected","Maternally immune"), lty=c(1,1),col=c(2,1), bty="n")

	plot(age.classes/12,(nt.store[m.inds,Tmax]+nt.store[s.inds,Tmax]+nt.store[i.inds1,Tmax]+nt.store[v.inds,Tmax])/diff(c(0,age.classes)), type="l", xlab="Ages", ylab="Tot pop over age")

    }
    return(list(Rt=Rt,nt.store=nt.store, m.inds=m.inds, s.inds= s.inds,i.indsv=i.indsv,s.indsv=s.indsv,i.inds1=i.inds1,s.inds1=s.inds1,i.inds2=i.inds2,s.inds2=s.inds2,i.inds3=i.inds3,s.inds3=s.inds3,
		i.inds4=i.inds4,v.inds=v.inds,
                age.classes=age.classes, Fmat=Fmat,Tmat=Tmat,store.beta=store.beta,waifw=waifw,burnin=burnin,Tmax=Tmax))

}


## function to work with the different age bin widths in terms of plotting
rescale <- function(x,ages) { 
	rc <- x/diff(c(0,ages))
	rc <- sum(x)*rc/sum(rc)
	return(rc) 
}


## rescale the Waifw to get a desired transmission magnitude
scale.Waifw <- function(R0, DFE, waifw){

    #more correct
    next.gen <- DFE*(1-exp(-waifw))

    #get the first eigen value
    cur.R0 <- Re(eigen(next.gen)$value[1])

    #More correct transform
    R.ratio <- R0/cur.R0; 
    waifw <- -log(1-R.ratio*(1-exp(-waifw)))

	#print(R.ratio)

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

get.smooth.WAIFW<-function(age.class.boundries = (1:120/12),
                           mu=12.71,sig=0.69, gam=0.17, delta=0){

    n.age.cats <- length(age.class.boundries)
    ages.to.use <- (age.class.boundries +
                    c(0,age.class.boundries[2:n.age.cats-1]))/2

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

get.polymod.WAIFW <- function (age.class.boundries = (1:90),
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
    n.age.cats <- length(age.class.boundries)
    lowpoints <-  c(0,age.class.boundries[2:n.age.cats-1])
    index <- (findInterval(lowpoints,ages.polymod.smooth));
    #set very large ages to all be the same as the largest age
    index[index>=length(ages.polymod.smooth)]=length(ages.polymod.smooth)-1

    #extract appropriate matrix, taking smoothed estimates for the ranges
    foi.matrix <- est$fhat[index+1,index+1]

    #adjust for fact that the width of your age class should not affect the number of contacts you make
    #foi.matrix <- foi.matrix/diff(c(0,age.class.boundries))

    colnames(foi.matrix) <- age.class.boundries
    rownames(foi.matrix) <- age.class.boundries
    return(foi.matrix)
}








#### Do some messing around  ############################################################################################################
exploreResults <- function(){

	## tweak this to adjust the rate of increase to get something close to equilbirium (corresponds to lambda=1)
	xx <- findStableStruct(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
				mort=c(rep(0,72),1),fert =  c(rep(0,66),rep(0.0007,7)))

	plot(xx$age.classes/12,xx$stable.age,xlab="age", ylab="", type="l")
	xx$lambda
	## I have picked out the survival / fertility I got from this for the examples that follow 
	## easier to control/understand if pop isn't growing like crazy, or shrinking. 

	## check out the time-series - set no vaccination 
	## with a flat waifw, putting beta into the waifw matrix at all points
	tmp <- iterateSIR(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), mort=c(rep(0,72),1),	
                          waifw=matrix(0.0002,73,73),
			alpha=0.4,
                         waning.maternal.par=0.5,recovery.infection=c(1,1,1,1),
			vaccination=c(rep(0,68),rep(0.001,5)),	
			contribution.to.trans=c(1,0,0,0),					
			time.step=1/4,
                         fert= c(rep(0,66),rep(0.0007,7)),
			burnin=(52*20),Tmax=(52*30), do.plot=TRUE,start.pop=c(),gamma=0.97)

	

	## now with lower seasonality to show the annual (rather than biennial) dynamics
	tmp1 <- iterateSIR(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), mort=c(rep(0,72),1),	
                          waifw=matrix(0.0002,73,73),
			alpha=0.1,
                         waning.maternal.par=0.5,recovery.infection=c(1,1,1,1),
			vaccination=c(rep(0,68),rep(0.001,5)),						
			contribution.to.trans=c(1,0,0,0),					
			time.step=1/4,
                         fert= c(rep(0,66),rep(0.0007,7)),
			burnin=(52*20),Tmax=(52*30), do.plot=TRUE,start.pop=c(),gamma=0.97)


	par(mfrow=c(1,2), bty="l")
  	plot((tmp$burnin:tmp$Tmax)/52, colSums(tmp$nt.store[tmp1$i.inds1,tmp1$burnin:tmp$Tmax]), type="l",col=2,xlab="Time", ylab="Total cases")
	points((tmp$burnin:tmp$Tmax)/52, colSums(tmp1$nt.store[tmp$i.inds1,tmp$burnin:tmp$Tmax]), type="l",col=1)
	abline(v=c(0:100),lty=1,col="grey")

	plot(tmp1$age.classes,tmp1$nt.store[tmp$i.inds1,52*30]/sum(tmp1$nt.store[tmp$i.inds1,52*30]), type="l", xlim=c(0,60), xlab="Age (months)", ylab="Relative propn primarny infected")
	points(tmp$age.classes,tmp$nt.store[tmp$i.inds1,52*30]/sum(tmp$nt.store[tmp$i.inds1,52*30]),type="l",col=2)
	legend("topright",legend=c("annual dynamics", "beinnial dynamics"),lty=1,col=c(1,2), bty="n")


 	plot((tmp$burnin:tmp$Tmax)/52, colSums(tmp$nt.store[tmp1$i.inds4,tmp1$burnin:tmp$Tmax]), type="l",col=2,xlab="Time", ylab="Total cases")
	points((tmp$burnin:tmp$Tmax)/52, colSums(tmp$nt.store[tmp$i.inds3,tmp$burnin:tmp$Tmax]), type="l",col="orange")
	points((tmp$burnin:tmp$Tmax)/52, colSums(tmp$nt.store[tmp$i.inds2,tmp$burnin:tmp$Tmax]), type="l",col=1)
	points((tmp$burnin:tmp$Tmax)/52, colSums(tmp$nt.store[tmp$i.inds1,tmp$burnin:tmp$Tmax]), type="l",col="grey")
	abline(v=c(0:100),lty=3,col="grey")




}


exploreReal <- function(){

	fname <- ""

	#run with Rachel's seasonality and pop structure  
	## download from: https://github.com/rebaker64/NPIs
	dfcases <- read.csv(paste(fname,"data/RSV_data_US.csv",sep=""))
	#make a timevariable that starts in year 0 just because that is how it is plotted
	times.cases <- dfcases$year+dfcases$week/52-2020+4

	load(paste(fname,"data/prelimRSVstatefitFlorida.RData",sep=""))
	print(mean(outall$sea_beta))

	# figure out what fertility should be to get these weekly births 
	chk <- findFert(target.births=mean(outall$births),fert.multiplier=seq(0.00001,0.02,length=10),popsize=outall$pop[1],
			age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
			mort=c(rep(1e-9,72),1), 
                         fert =  c(rep(0,66),rep(0.1,7)), time.step=1/4)

	tmp <- iterateSIR(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), mort=c(rep(1e-9,72),1),	
                          waifw=matrix(median(outall$sea_beta),73,73),
			alpha=0.2,#note that this does nothing in this formulation
			discrete.seas=outall$sea_beta,
                         waning.maternal.par=0.5,
			contribution.to.trans=c(1,0,0,0),					
			vaccination=c(rep(0,68),rep(0.001,5)),						
			time.step=1/4,
                         fert=  c(rep(0,66),rep(0.1,7))*chk,
			burnin=(52*20),Tmax=(52*30), do.plot=TRUE,start.pop=c(),start.pop.size=outall$pop[1],gamma=0.97)

	load(paste(fname,"data/prelimRSVstatefitTEXAS.RData",sep=""))
	print(mean(outall$sea_beta))
	# figure out what fertility should be to get these weekly births 
	chk <- findFert(target.births=mean(outall$births),fert.multiplier=seq(0.00001,0.02,length=10),popsize=outall$pop[1],
			age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
			mort=c(rep(1e-9,72),1), 
                         fert =  c(rep(0,66),rep(0.1,7)), time.step=1/4)

	tmp1 <- iterateSIR(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), mort=c(rep(1e-9,72),1),	
                         waifw=matrix(median(outall$sea_beta),73,73),
			alpha=0.2,
			discrete.seas=outall$sea_beta,
                         waning.maternal.par=0.5,
			contribution.to.trans=c(1,0,0,0),					
			vaccination=c(rep(0,68),rep(0.001,5)),						
			time.step=1/4,
                         fert=  c(rep(0,66),rep(0.1,7))*chk,
			burnin=(52*20),Tmax=(52*30), do.plot=TRUE,start.pop=c(),start.pop.size=outall$pop[1],gamma=0.97)


	par(mfrow=c(1,3))
	chs <- dfcases$state=="Florida"
	plot(times.cases[chs],(dfcases$percent_specimen_positive[chs])/max(dfcases$percent_specimen_positive[chs]),type="b", xlab="", ylab="cases", lty=1,lwd=2,pch=19)
 	points(((tmp$burnin:tmp$Tmax)-tmp$burnin)/52, colSums(tmp$nt.store[tmp$i.inds1,tmp$burnin:tmp$Tmax])/max(colSums(tmp$nt.store[tmp$i.inds1,tmp$burnin:tmp$Tmax])), type="l", col="blue")
	abline(v=c(0:100),lty=3)
	title("Florida")
		
	chs <- dfcases$state=="Texas"
	plot(times.cases[chs],(dfcases$percent_specimen_positive[chs])/max(dfcases$percent_specimen_positive[chs]),type="b", xlab="", ylab="cases", lty=1,lwd=2,pch=19)
	points(((tmp$burnin:tmp$Tmax)-tmp$burnin)/52, colSums(tmp1$nt.store[tmp1$i.inds1,tmp$burnin:tmp$Tmax])/max(colSums(tmp1$nt.store[tmp$i.inds1,tmp$burnin:tmp$Tmax])), 
		type="l", col="blue")
	abline(v=c(0:100),lty=3)
	title("Texas")
	
	#plot(tmp1$age.classes,tmp1$nt.store[tmp$i.inds1,52*30]/sum(tmp1$nt.store[tmp$i.inds1,52*30]), type="l", xlim=c(0,60),xlab="Age (months)", ylab="Relative propn infected", ylim=c(0,0.01))
	#points(tmp$age.classes,tmp$nt.store[tmp$i.inds1,52*30]/sum(tmp$nt.store[tmp$i.inds1,52*30]),type="l",col=2)

	indexes <- seq((52*29),(52*30),by=4)
	cols1 <- colorRampPalette(c(rgb(0,0,1,1), rgb(0,0,1,0)), alpha = TRUE)(7); cols1 <- c(cols1,cols1[6:1])
	cols2 <- colorRampPalette(c(rgb(1,0,0,1), rgb(1,0,0,0)), alpha = TRUE)(7); cols2 <- c(cols2,cols2[6:1])

	matplot(tmp1$age.classes,tmp1$nt.store[tmp$i.inds1,indexes]/sum(tmp1$nt.store[tmp$i.ind1s,52*30]), type="l",lty=1,
		xlim=c(0,58), xlab="Age (months)", ylab="Relative propn infected", ylim=c(0,0.02),col=cols2)
	matplot(tmp$age.classes,tmp$nt.store[tmp$i.inds1,indexes]/sum(tmp$nt.store[tmp$i.inds1,52*30]),type="l",lty=1,col=cols1, add=TRUE)
	legend("topright",legend=c("Florida", "Texas"),lty=1,col=c(4,2), bty="n")

	points(tmp$age.classes,tmp1$nt.store[tmp$i.inds1,(52-16)*30]/sum(tmp1$nt.store[tmp$i.inds1,(52-16)*30]),type="l",col=1,lwd=2,lty=3)
	points(tmp$age.classes,tmp$nt.store[tmp$i.inds1,(52-16)*30]/sum(tmp$nt.store[tmp$i.inds1,(52-16)*30]),type="l",col=1,lwd=2,lty=3)


}




exploreFunkyWAIFW <- function(){

	fname <- ""
	load(paste(fname,"data/prelimRSVstatefitFlorida.RData",sep=""))

	#pick a waifw and glance at it (if the transmission is all way too low or too old, the R0 scaling won't work) 
	##theoretical 
	waifw.mat <- get.smooth.WAIFW(age.class.boundries=c(1:60,seq(72,120,by=12),seq(180,600,by=60))/12)
	image(c(1:60,seq(72,120,by=12),seq(180,600,by=60))/12,c(1:60,seq(72,120,by=12),seq(180,600,by=60))/12,waifw.mat,xlab="age (years)", ylab="age contact (years)")
	##empirical (polymod)
	waifw.mat <- get.polymod.WAIFW(age.class.boundries = c(1:60,seq(72,120,by=12),seq(180,600,by=60))/12,
                               country="GB",bandwidth=c(3,3),do.touch=FALSE) 
	image(c(1:60,seq(72,120,by=12),seq(180,600,by=60))/12,c(1:60,seq(72,120,by=12),seq(180,600,by=60))/12,waifw.mat,xlab="age (years)", ylab="age contact (years)")
	

	#get stable age structure
	xx <- findStableStruct(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
				mort= c(rep(0,72),1),fert =  c(rep(0,66),rep(1,7)))

	# figure out what fertility should be to get these weekly births 
	chk <- findFert(target.births=mean(outall$births),fert.multiplier=seq(0.00001,0.002,length=10),popsize=outall$pop[1],
			age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
			mort=c(rep(1e-9,72),1), 
                         fert =  c(rep(0,66),rep(1,7)), time.step=1/4)

	#get stable age structure with chk in place
	xx <- findStableStruct(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
				mort= c(rep(0,72),1),fert =  c(rep(0,66),rep(1,7))*chk)

	## rescale the waifw
	scaled.waifw.mat <- scale.Waifw(R0=median(outall$sea_beta)*outall$pop[1], DFE=xx$stable.age*outall$pop[1], waifw=waifw.mat)

	tmp <- iterateSIR(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), mort=c(rep(0,72),1),	
                         waifw=scaled.waifw.mat,
			alpha=0.2,
			discrete.seas=outall$sea_beta,
                         waning.maternal.par=0.5,
			contribution.to.trans=c(1,0,0,0),					
			vaccination=c(rep(0,68),rep(0.001,5)),						
			time.step=1/4,
                         fert=  c(rep(0,66),rep(1,7))*chk,
			burnin=(52*20),Tmax=(52*30), do.plot=TRUE,start.pop=c(),start.pop.size=outall$pop[1],gamma=0.97)


	dfcases <- read.csv(paste(fname,"data/RSV_data_US.csv",sep=""))
	times.cases <- dfcases$year+dfcases$week/52-2020+4

	par(mfrow=c(1,3))
	chs <- dfcases$state=="Florida"
	plot(times.cases[chs],(dfcases$percent_specimen_positive[chs])/max(dfcases$percent_specimen_positive[chs]),type="b", xlab="", ylab="cases", lty=1,lwd=2,pch=19)
 	points(((tmp$burnin:tmp$Tmax)-tmp$burnin)/52, colSums(tmp$nt.store[tmp$i.inds1,tmp$burnin:tmp$Tmax])/max(colSums(tmp$nt.store[tmp$i.inds1,tmp$burnin:tmp$Tmax])), 
		type="l", col="blue")
	abline(v=c(0:100),lty=3)
	title("Florida")



}


