
#setwd("/Users/cmetcalf/Dropbox (Princeton)/RSV-Age-structure/")

## Function to bring in this chunk of code. 
#
reload.source <- function(fname=""){
    source(paste(fname,"source/AgeStructModel-RSV.R",sep=""))
}



## Function to build survival, aging and infectious transitions (T matrix called 'U' in the methods)
#
buildTMatrixSlowerDemog <- function(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)),  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years (could extend)
                                    mort=c(rep(1e-9,72),1),						#neglgible death till age 50 (could make more reaslitic)
                                    beta=c(0.1),							#transmission (can be single value or vector of same length as age)
                                    waning.maternal.par=0.5,					#rate of waning of maternal immunity
                                    waning.immunity=0, 						#rate of waning of general immunity
				                            vaccination=c(rep(0,68),rep(0.5,5)),				#vector of length ages (e.g., vaccinate at age 20, vaccinate 50% of the population)		
				                            time.step=1/4){						#time-step is 1/4 of a month, aka 1 week (note slight weirdness with month/year 12 months = 52 weeks)
 
    ## matrix dimensions, aging rate 
    nage <- length(age.classes)
    aging.rate <- time.step/diff(c(0,age.classes))

    ## adjustments for cases where ebta or waning immunity are constant across age and only one number provided
    if (length(beta)==1) beta <- rep(beta,nage)
    if (length(waning.immunity)==1) waning.immunity <- rep(waning.immunity,nage)

    ## waning of maternal immunity: first calculate the overall pattern, i.e., proportion of individuals having 'survived' waning in every age class	
    pmaternal <- exp(-c(0,age.classes)*waning.maternal.par)
    ## What leaks out specifically in each age class - taking bottom and top of age class and getting proportion that have gone
    waning.maternal <- (pmaternal[1:(length(pmaternal)-1)]-pmaternal[2:length(pmaternal)])/pmaternal[1:(length(pmaternal)-1)]	

    #two matrices for storage	    
    mat1 <- matrix(0,4,4)
    Tmat <- matrix(0,4*nage,4*nage)

    #loop over ages
    for (j in 1:nage) {
        #fill in epi matrix	
        mat1[] <- 0
        mat1[1,1] <- 1-waning.maternal[j]
        mat1[2,1] <- waning.maternal[j]
        mat1[2,2] <- (1-beta[j])*(1-vaccination[j])
        mat1[3,2] <- beta[j]*(1-vaccination[j])
        mat1[4,2] <- vaccination[j]  	## assuming that vaccination precedes transmission 
        mat1[4,3] <- 1
        mat1[4,4] <- 1-waning.immunity[j]
        mat1[2,4] <- waning.immunity[j]

        #put in surv
        surv <- rep(1-mort[j],4);
   
        #fill in Tmatrix
        if (j!=nage) {
            #go into next age class
            Tmat[(j*4+1):(j*4+4),((j-1)*4+1):(j*4)] <- mat1*surv*aging.rate[j]

            #stay in this age class - only lose maternal immunity when change an age class
            mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
            Tmat[((j-1)*4+1):(j*4),((j-1)*4+1):(j*4)] <- mat2*surv*(1-aging.rate[j])

        } else {
            #stay in the last age class
            Tmat[((j-1)*4+1):(j*4),((j-1)*4+1):(j*4)] <- mat1*surv
        }}

	return(Tmat)

}


## Funtion to build fertility transitions
#
buildFMatrix <- function(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)),  
                         fert =  c(rep(0,66),rep(0.1,7))){ 	#one for every age class. start reproducing at age 20, and here assume ~ flat
	nage <- length(age.classes)

	Fmat <- matrix(0,4*nage,4*nage)
	for (j in 1:nage) {
		Fmat[1,((j-1)*4+1):(j*4)] <- c(0,0,fert[j],fert[j]) 	#the maternally immune born from infected and recovered
		Fmat[2,((j-1)*4+1):(j*4)]<- c(fert[j],fert[j],0,0)  	#the susceptible (mom didn't get sick, kid born susceptible)
	}

	return(Fmat)
}


## Function to find stable age structure; population rate of increase, etc; from fertility and mortality profiles (no classification by epidemiological class)

findStableStruct <- function(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
			                      mort=c(rep(1e-9,72),1), fert =  c(rep(0,66),rep(0.1,7)), time.step=1/4){


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

	return(list(stable.age=stable.age,lambda=lambda,reprod.value=reprod.value,age.classes=age.classes,Fmat=Fmat,Tmat=Tmat))
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
                        mort=c(rep(0,72),1),					        #neglgible death till age 50 (could make more reaslitic)
                        waifw=matrix(0.0002,73,73),					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
			alpha=0.2,seas.offset=0,							#parameter governing sine wave seasonality in transmission	
			discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                        waning.maternal.par=0.5,						#rate of waning of maternal immunity
                        waning.immunity=0, 						#rate of waning of immunity in general
			vaccination=c(rep(0,68),rep(0.5,5)),				#vaccination strategy (e.g., this default: vaccinate at age 20, vaccinate 50% of the population)		
		        time.step=1/4,							#time-step (1/4 of a month = 1 week)
                        fert= c(rep(0,66),rep(0.0008,7)),					#pattern of fertility over age
			alpha.births = 0,						#pattern of seasonality in births	
			burnin=10,Tmax=150, 
			do.plot=FALSE,	
		        start.pop=c(),							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
			start.pop.size=100000,						#starting total population size
			start.full.pop.struct=c(),					#starting full pop structure - for if you want to switch off vacc
			vacc.t.multip=52*10,						#time to multiply out the vacc by chosen value
			value.multip.vacc=1,						#value for 
			value.add.vacc=0,
			v.times=c(),              					#if you implement this, the vaccination vector above should be 1 in age classes where vaccination occurs
			effect.npis=1,							#function to multiply beta by to rflect NPIs, should be of same length as Tmax if in use
			gamma=0.97, do.unity=FALSE){					#this parameter mitigates explosiveness emerging from discrete time.
 

    #indexing - the rows for maternal, susceptible, etc
    m.inds <- seq(1,length(age.classes)*4,by=4)
    s.inds <- m.inds+1
    i.inds <- m.inds+2
    r.inds <- m.inds+3

    #get Fmatrix (basic won't change - might multiply by seasonal flucutation, but can happen as single parameter)
    Fmat <- buildFMatrix(age.classes=age.classes,fert=fert)

    
    #over-write everything if have a starting point 
    if (length(start.full.pop.struct)!=0) { 
      nt <- start.full.pop.struct
    } else {
    #initiate populations
    find.expected.stable <- findStableStruct(age.classes=age.classes,mort= mort,fert =  fert, time.step=time.step)
    if (length(start.pop)==0) start.pop <- start.pop.size*find.expected.stable$stable.age

    nt<-rep(0,length(age.classes)*4)
    nt[s.inds] <- start.pop  			 #susceptibles
    nt[i.inds[1:min(60,length(age.classes))]] <- 10 			 #seed some infected (otherwise infection can't take off)
    nt[s.inds[1:min(60,length(age.classes))]] <- nt[s.inds[1:min(60,length(age.classes))]]-10 	 #remove to keep pop size the same
    }

    #print(nt)

    ## storage    
    Rt<-store.beta <- store.beta.unity <- store.births <- rep(NA,Tmax)		 #storage
    nt.store <- matrix(NA,length(nt),Tmax)

    ## modulate vaccination
    if (length(v.times)==0) { 
      v.times <- rep(value.multip.vacc,Tmax)	#fill with the multiplier 
      v.times[1:vacc.t.multip] <- 1         	#make it 1 before the switch point
    }
    v.add <- rep(0,Tmax) 
    v.add[vacc.t.multip:Tmax] <- value.add.vacc


    ## modulate NPI effects
    if (length(effect.npis)==1) effect.npis <- rep(effect.npis,Tmax)
    if (length(effect.npis)<Tmax) effect.npis <- c(effect.npis,rep(1,Tmax-length(effect.npis)))
  
    
    #initiate pop
    nt.store[,1] <-  nt

    #seasonal transmission - each has a mean / median of 1 TODO this is not general (26, 52)
    if (length(discrete.seas)==0) seas <- pmax(1+alpha*cos(2*pi*(seas.offset+(1:26))/26),0)	else { seas <- discrete.seas/median(discrete.seas); beta <-1 }

    start.pop.size <- sum(nt) #this for where it isn't defined before and you need it to constrain pop size from getting huge
    
    #TODO quick fix for adding vaccination; assume young age classes in months and vaccination at 12 months 
    age.struct.vacc <- rep(0,length(age.classes))
    age.struct.vacc[12] <- 1

    for (t in 1:Tmax) {

      	#extract rate of transmission to everage age TODO
	      phi <- (1-exp(-effect.npis[t]*seas[1+t%%26]*waifw%*%(nt[i.inds]^gamma))); 
	      beta.trans <- t(phi); 

        #build the T matrix 
        Tmat <- buildTMatrixSlowerDemog(age.classes=age.classes,
                                        mort=mort,   	
                                        beta=beta.trans,
                                        waning.maternal.par=waning.maternal.par,
                                        waning.immunity=waning.immunity,
                                        vaccination=(vaccination*v.times[t]+age.struct.vacc*v.add[t]), time.step=time.step)

        #print(waning.maternal.par)
        #print(Tmat[1:10,1:10])
        #print(which(is.na(Tmat),arr.ind=TRUE))[1:15,]
        
	      #introduce seasonality in births (if alpha.births=0, this will just be x1)
	      seas.births <- (1+alpha.births*cos(2*pi*t/26))	

	      #use matrix multiplication to get the population structure at the next time-step
	      nt1<-(Tmat+Fmat*seas.births) %*% nt
	      #print(c(t,range(Tmat), sum(nt1)))
	      

	##Introduce effective transmission [see Mahmud thesis]
	if (do.unity) store.beta.unity[t] <- (1/((sum(nt[i.inds]^gamma))*sum(nt[s.inds])))*sum(beta.trans*nt[s.inds])
	if (do.unity) store.births[t] <- sum((Fmat*seas.births) %*% nt)

        sum.nt1<-sum(nt1)
        Rt[t]<-log(sum.nt1)	#pop growth rate
        nt<-nt1 			      #update the population 
        nt.store[,t] <- nt	#store pop
        store.beta[t] <- max(beta.trans)	#store transmission


        #if earlier than burnin, and no infecteds, introduce some
        nt[i.inds] <- pmin(nt[i.inds],start.pop.size)
       # if (t>300) print(c(range(phi),sum(nt[i.inds])))
       # if (t>300) print(seas[1+t%%26])
        if (sum(nt[i.inds])<1 & t<burnin) {nt[i.inds] <- 1/length(age.classes); nt[s.inds] <- pmax(nt[s.inds]-1/length(age.classes),0) }
    }
 
    if (do.plot) { ## mostly sanity checks / diagnoses
        par(mfrow=c(2,1))
        image((burnin:Tmax)/52,age.classes/12,log(t(nt.store[s.inds,burnin:Tmax])+1), xlab="Time", ylab="Age", main="susceptible", ylim=c(0,15))
        contour((burnin:Tmax)/52,age.classes/12,t(nt.store[s.inds,burnin:Tmax])+1,add=TRUE)
        image((burnin:Tmax)/52,age.classes/12,log(t(nt.store[i.inds,burnin:Tmax])+1), xlab="Time", ylab="Age", main="infected", ylim=c(0,15))
        contour((burnin:Tmax)/52,age.classes/12,t(nt.store[i.inds,burnin:Tmax])+1,add=TRUE)


	par(mfrow=c(2,2))
  	plot((burnin:Tmax)/52, colSums(nt.store[i.inds,burnin:Tmax]), type="l",col=1,xlab="Time", ylab="Total cases")
	abline(v=c(0:100),lty=3)
	plot((burnin:Tmax)/52, colSums(nt.store[m.inds,burnin:Tmax])+colSums(nt.store[s.inds,burnin:Tmax])+
				colSums(nt.store[i.inds,burnin:Tmax])+colSums(nt.store[r.inds,burnin:Tmax]), type="l",col=1,xlab="Time", ylab="Total pop size")#, xlim=c(120,134))
	
	plot(age.classes/12,nt.store[m.inds,Tmax]/diff(c(0,age.classes)), type="l", xlab="Ages", ylab="Number", xlim=c(0,10))
	points(age.classes/12,nt.store[m.inds,burnin]/diff(c(0,age.classes)), type="l")

	points(age.classes/12,nt.store[i.inds,Tmax]/diff(c(0,age.classes)),type="l",col=2)
	points(age.classes/12,nt.store[i.inds,burnin]/diff(c(0,age.classes)),col=2, type="l")
	legend("topright",legend=c("Infected","Maternally immune"), lty=c(1,1),col=c(2,1), bty="n")

	plot(age.classes/12,(nt.store[m.inds,Tmax]+nt.store[s.inds,Tmax]+nt.store[i.inds,Tmax]+nt.store[r.inds,Tmax])/diff(c(0,age.classes)), type="l", xlab="Ages", ylab="Tot pop over age")
	#check equilibrium -
	#plot(age.classes/12,(start.pop)/diff(c(0,age.classes)), type="l", xlab="Ages", ylab="Number")
	#points(age.classes/12,(nt.store[m.inds,1]+nt.store[s.inds,1]+nt.store[i.inds,1]+nt.store[r.inds,1])/diff(c(0,age.classes)), type="l",col=2)
	#points(age.classes/12,(nt.store[m.inds,10]+nt.store[s.inds,10]+nt.store[i.inds,10]+nt.store[r.inds,10])/diff(c(0,age.classes)), type="l",col=4)
	#points(age.classes/12,(nt.store[m.inds,50]+nt.store[s.inds,50]+nt.store[i.inds,50]+nt.store[r.inds,50])/diff(c(0,age.classes)), type="l",col=3)

    }
    return(list(Rt=Rt,nt.store=nt.store, i.inds=i.inds,s.inds=s.inds,m.inds=m.inds,r.inds=r.inds,
                age.classes=age.classes, Fmat=Fmat,Tmat=Tmat,store.beta=store.beta,store.beta.unity=store.beta.unity, store.births= store.births,
		waifw=waifw,burnin=burnin,Tmax=Tmax,v.times=v.times,vacc.t.multip=vacc.t.multip))

}







# Function to iterate out to equilibrium - but this is specifcally designed to mess with honemoons 
#
#
iterateSIRvacc <- function(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)),  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                        mort=c(rep(0,72),1),					        #neglgible death till age 50 (could make more reaslitic)
                        waifw=matrix(0.2,73,73),					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
			alpha=0.2,seas.offset=0,							#parameter governing sine wave seasonality in transmission	
			discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                        waning.maternal.par=0.5,						#rate of waning of maternal immunity
                        waning.immunity=0, 						#rate of waning of immunity in general
                        vacc.cover.start=0.9, vacc.cover.release=0.4,
		        time.step=1/4,							#time-step (1/4 of a month = 1 week)
                        fert= c(rep(0,72),1),					#pattern of fertility over age
			alpha.births = 0,						#pattern of seasonality in births	
			burnin=10,Tmax=(26*500), 
			do.plot=FALSE,	
		        start.pop=c(),							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
			start.pop.size=100000,						#starting total population size
			start.full.pop.struct=c(),					#starting full pop structure - for if you want to switch off vacc
			effect.npis=1,							#function to multiply beta by to rflect NPIs, should be of same length as Tmax if in use
			gamma=0.97, do.unity=FALSE){					#this parameter mitigates explosiveness emerging from discrete time.
 

    #indexing - the rows for maternal, susceptible, etc
    m.inds <- seq(1,length(age.classes)*4,by=4)
    s.inds <- m.inds+1
    i.inds <- m.inds+2
    r.inds <- m.inds+3

    #get Fmatrix (basic won't change - might multiply by seasonal flucutation, but can happen as single parameter)
    Fmat <- buildFMatrix(age.classes=age.classes,fert=fert)

    
    #over-write everything if have a starting point 
    if (length(start.full.pop.struct)!=0) { 
      nt <- start.full.pop.struct
    } else {
    #initiate populations
    find.expected.stable <- findStableStruct(age.classes=age.classes,mort= mort,fert =  fert, time.step=time.step)
    if (length(start.pop)==0) start.pop <- start.pop.size*find.expected.stable$stable.age

    nt<-rep(0,length(age.classes)*4)
    nt[s.inds] <- start.pop  			 #susceptibles
    nt[i.inds[1:min(60,length(age.classes))]] <- 10 			 #seed some infected (otherwise infection can't take off)
    nt[s.inds[1:min(60,length(age.classes))]] <- nt[s.inds[1:min(60,length(age.classes))]]-10 	 #remove to keep pop size the same
    }

    #print(nt)

    ## storage    
    Rt<-store.beta <- store.beta.unity <- store.births <- rep(NA,Tmax)		 #storage
    nt.store <- matrix(NA,length(nt),Tmax)

    ## modulate NPI effects
    if (length(effect.npis)==1) effect.npis <- rep(effect.npis,Tmax)
    if (length(effect.npis)<Tmax) effect.npis <- c(effect.npis,rep(1,Tmax-length(effect.npis)))
  
    ## set up vaccination    
    vacc.vals <- rep(0,Tmax)
    vacc.vals[1:(85*26+3)] <- 0
    vacc.vals[(85*26+3+1):Tmax] <- vacc.cover.start
    vacc.vals[(85*26+3+1+26*5+1):Tmax] <- vacc.cover.release

    #initiate pop
    nt.store[,1] <-  nt

    #seasonal transmission - each has a mean / median of 1 TODO this is not general (26, 52)
    if (length(discrete.seas)==0) seas <- pmax(1+alpha*cos(2*pi*(seas.offset+(1:26))/26),0)	else { seas <- discrete.seas/median(discrete.seas); beta <-1 }

    start.pop.size <- sum(nt) #this for where it isn't defined before and you need it to constrain pop size from getting huge
    
    #TODO quick fix for adding vaccination; assume young age classes in months and vaccination at 12 months 
    age.struct.vacc <- rep(0,length(age.classes))
    age.struct.vacc[12:15] <- 1 #note that previous just has till 13

    for (t in 1:Tmax) {

      	#extract rate of transmission to everage age TODO
	      phi <- (1-exp(-effect.npis[t]*seas[1+t%%26]*waifw%*%(nt[i.inds]^gamma))); 
	      beta.trans <- t(phi); 

        #build the T matrix 
        Tmat <- buildTMatrixSlowerDemog(age.classes=age.classes,
                                        mort=mort,   	
                                        beta=beta.trans,
                                        waning.maternal.par=waning.maternal.par,
                                        waning.immunity=waning.immunity,
                                        vaccination=vacc.vals[t]*age.struct.vacc, time.step=time.step)

        #print(waning.maternal.par)
        #print(Tmat[1:10,1:10])
        #print(which(is.na(Tmat),arr.ind=TRUE))[1:15,]
        
	      #introduce seasonality in births (if alpha.births=0, this will just be x1)
	      seas.births <- (1+alpha.births*cos(2*pi*t/26))	

	      #use matrix multiplication to get the population structure at the next time-step
	      nt1<-(Tmat+Fmat*seas.births) %*% nt
	      #print(c(t,range(Tmat), sum(nt1)))
	      

	##Introduce effective transmission [see Mahmud thesis]
	if (do.unity) store.beta.unity[t] <- (1/((sum(nt[i.inds]^gamma))*sum(nt[s.inds])))*sum(beta.trans*nt[s.inds])
	if (do.unity) store.births[t] <- sum((Fmat*seas.births) %*% nt)

        sum.nt1<-sum(nt1)
        Rt[t]<-log(sum.nt1)	#pop growth rate
        nt<-nt1 			      #update the population 
        nt.store[,t] <- nt	#store pop
        store.beta[t] <- max(beta.trans)	#store transmission


        #if earlier than burnin, and no infecteds, introduce some
        nt[i.inds] <- pmin(nt[i.inds],start.pop.size)
       # if (t>300) print(c(range(phi),sum(nt[i.inds])))
       # if (t>300) print(seas[1+t%%26])
        if (sum(nt[i.inds])<1 & t<burnin) {nt[i.inds] <- 1/length(age.classes); nt[s.inds] <- pmax(nt[s.inds]-1/length(age.classes),0) }
    }
 
    if (do.plot) { ## mostly sanity checks / diagnoses
        par(mfrow=c(2,1))
        image((burnin:Tmax)/52,age.classes/12,log(t(nt.store[s.inds,burnin:Tmax])+1), xlab="Time", ylab="Age", main="susceptible", ylim=c(0,15))
        contour((burnin:Tmax)/52,age.classes/12,t(nt.store[s.inds,burnin:Tmax])+1,add=TRUE)
        image((burnin:Tmax)/52,age.classes/12,log(t(nt.store[i.inds,burnin:Tmax])+1), xlab="Time", ylab="Age", main="infected", ylim=c(0,15))
        contour((burnin:Tmax)/52,age.classes/12,t(nt.store[i.inds,burnin:Tmax])+1,add=TRUE)


	par(mfrow=c(2,2))
  	plot((burnin:Tmax)/52, colSums(nt.store[i.inds,burnin:Tmax]), type="l",col=1,xlab="Time", ylab="Total cases")
	abline(v=c(0:100),lty=3)
	plot((burnin:Tmax)/52, colSums(nt.store[m.inds,burnin:Tmax])+colSums(nt.store[s.inds,burnin:Tmax])+
				colSums(nt.store[i.inds,burnin:Tmax])+colSums(nt.store[r.inds,burnin:Tmax]), type="l",col=1,xlab="Time", ylab="Total pop size")#, xlim=c(120,134))
	
	plot(age.classes/12,nt.store[m.inds,Tmax]/diff(c(0,age.classes)), type="l", xlab="Ages", ylab="Number", xlim=c(0,10))
	points(age.classes/12,nt.store[m.inds,burnin]/diff(c(0,age.classes)), type="l")

	points(age.classes/12,nt.store[i.inds,Tmax]/diff(c(0,age.classes)),type="l",col=2)
	points(age.classes/12,nt.store[i.inds,burnin]/diff(c(0,age.classes)),col=2, type="l")
	legend("topright",legend=c("Infected","Maternally immune"), lty=c(1,1),col=c(2,1), bty="n")

	plot(age.classes/12,(nt.store[m.inds,Tmax]+nt.store[s.inds,Tmax]+nt.store[i.inds,Tmax]+nt.store[r.inds,Tmax])/diff(c(0,age.classes)), type="l", xlab="Ages", ylab="Tot pop over age")
	#check equilibrium -
	#plot(age.classes/12,(start.pop)/diff(c(0,age.classes)), type="l", xlab="Ages", ylab="Number")
	#points(age.classes/12,(nt.store[m.inds,1]+nt.store[s.inds,1]+nt.store[i.inds,1]+nt.store[r.inds,1])/diff(c(0,age.classes)), type="l",col=2)
	#points(age.classes/12,(nt.store[m.inds,10]+nt.store[s.inds,10]+nt.store[i.inds,10]+nt.store[r.inds,10])/diff(c(0,age.classes)), type="l",col=4)
	#points(age.classes/12,(nt.store[m.inds,50]+nt.store[s.inds,50]+nt.store[i.inds,50]+nt.store[r.inds,50])/diff(c(0,age.classes)), type="l",col=3)

    }
    return(list(Rt=Rt,nt.store=nt.store, i.inds=i.inds,s.inds=s.inds,m.inds=m.inds,r.inds=r.inds,
                age.classes=age.classes, Fmat=Fmat,Tmat=Tmat,store.beta=store.beta,store.beta.unity=store.beta.unity, store.births= store.births,
		waifw=waifw,burnin=burnin,Tmax=Tmax, vacc.vals = vacc.vals,seas=seas))

}




## preliminary function to explore vaccination coverage effect on average age of infection
exploreVaccination <- function(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)),  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                          mort=c(rep(0,72),1),						#neglgible death till age 50 (could make more reaslitic)
                     	waifw=matrix(0.01,73,73),
			alpha=0.2,							#seasonality in transmission	
                         waning.maternal.par=0.5,
                         waning.immunity=0, 
			time.step=1/4,
                         fert= c(rep(0,66),rep(0.0007,7)),
			vaccination.template=c(rep(0,66),rep(1,7)),		#vaccinate older individuals - hacky variant, only hash out opeions!
			#vaccination.template=rep(1,73),				#vaccinate kids
			#vaccination.template=c(0,1,rep(0,71)),			#vaccinate kids at 2 months old
			burnin=(52*15),Tmax=(52*30),start.pop=c(),gamma=0.9,n.test=5, do.plot=FALSE){


	vacc.test <- seq(0,0.95,length=n.test) 
	total <- scaled.store.recovered <- scaled.store.susceptible <- scaled.store.infected <- scaled.store.maternally.immune <- matrix(NA,length(vacc.test),length(age.classes))
	age.peak.infected <- avg.age.infection <- number.infected <- rep(NA,n.test)
	mid.age.classes <- c(c(0,age.classes[-length(age.classes)])+age.classes)*0.5
	#print(mid.age.classes)

	for (j in 1:length(vacc.test)) { 

		print(j)
		tmp <- iterateSIR(age.classes=age.classes,mort=mort,waifw=waifw,alpha=alpha,waning.maternal.par=waning.maternal.par,
			waning.immunity=waning.immunity,
			vaccination=vaccination.template*vacc.test[j],
			time.step=time.step,fert=fert,burnin=burnin,Tmax=Tmax,do.plot=FALSE,start.pop=start.pop,gamma=gamma)
		scaled.store.infected[j,] <- rescale(tmp$nt.store[tmp$i.inds,Tmax],tmp$age.classes)
		scaled.store.maternally.immune[j,] <- rescale(tmp$nt.store[tmp$m.inds,Tmax],tmp$age.classes)
		scaled.store.susceptible[j,] <- rescale(tmp$nt.store[tmp$s.inds,Tmax],tmp$age.classes)
		scaled.store.recovered[j,] <- rescale(tmp$nt.store[tmp$r.inds,Tmax],tmp$age.classes)
		total[j,] <- rescale(tmp$nt.store[tmp$m.inds,Tmax]+tmp$nt.store[tmp$s.inds,Tmax]+tmp$nt.store[tmp$i.inds,Tmax]+tmp$nt.store[tmp$r.inds,Tmax],tmp$age.classes)
				
	
		number.infected[j] <- sum(tmp$nt.store[tmp$i.inds,(Tmax-52):Tmax])
		chs <- age.classes<120 ## focus on averages and peaks in younger children
		avg.age.infection[j] <- sum(mid.age.classes[chs]*tmp$nt.store[tmp$i.inds[chs],Tmax]/sum(tmp$nt.store[tmp$i.inds[chs],Tmax]))/12
		age.peak.infected[j] <- mid.age.classes[scaled.store.infected[j,]==max(scaled.store.infected[j,])]/12

		
	}
	
	if (do.plot) { 

	par(mfrow=c(1,5))
	cols <- colorRampPalette(c("red", "blue"))(length(vacc.test))
	mid.point.age <- 0.5*(c(0,tmp$age.classes[-length(tmp$age.classes)])+tmp$age.classes)/12
	matplot(mid.point.age,t(scaled.store.maternally.immune), type="l",col=cols, xlab="Age (years)", ylab="Number maternally immune", lty=1, xlim=c(0,5))
	matplot(mid.point.age,t(scaled.store.susceptible), type="l",col=cols, xlab="Age (years)", ylab="Number susceptible", lty=1, xlim=c(0,5))
	matplot(mid.point.age,t(scaled.store.infected), type="l",col=cols, xlab="Age (years)", ylab="Number infected",  lty=1, xlim=c(0,5))
	matplot(mid.point.age,t(scaled.store.recovered), type="l",col=cols, xlab="Age (years)", ylab="Number recovered", lty=1)
	matplot(mid.point.age,t(total), type="l",col=cols, xlab="Age (years)", ylab="Total",  lty=1)

	}

	return(list(tmp=tmp, age.classes=age.classes,scaled.store.infected=scaled.store.infected,scaled.store.maternally.immune=scaled.store.maternally.immune,total=total, 
			avg.age.infection = avg.age.infection, number.infected=number.infected,age.peak.infected=age.peak.infected,vacc.test=vacc.test))

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
    #image(next.gen)
    #print(range(next.gen))
    #print(range(waifw))
    #print(range(DFE))
    #print("done")
    next.gen[is.na(next.gen)] <- mean(next.gen,na.rm=TRUE)
    next.gen[!is.finite(next.gen)] <- max(next.gen,na.rm=TRUE)
    
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
  browser()

   #bring in the polymod data on contacts for UK
    # polymod <- read.csv("data/polymodRaw.csv")
    polymod <- epimdr2::polymod
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





extractAge <- function(mat,ages=c(1:60,seq(72,120,by=12),seq(180,600,by=60))) { 

	mid.age <- 0.5*(c(0,ages[-length(ages)])+ages)
	avg.age <- var.age <- rep(NA,ncol(mat))	

	for (t in 1:ncol(mat)) { 
		pt <- mat[,t]/sum(mat[,t])
		#plot(mid.age/12,pt, type="l")
		avg.age[t] <- sum(pt*mid.age)
		var.age[t] <- sum(pt*(mid.age-avg.age[t])^2)		
	}

	return(list(avg.age=avg.age,var.age=var.age))

}

findRate <- function(par,ages,prop){
    pred <- 1-exp(-par*ages)
    u <- 1/sum((pred-prop)^2)
    return(u)
}

#update age using cumulative and then parametric form to get at average time
## this is not super different from the basic 
extractAge1 <- function(mat,ages=c(1:60,seq(72,120,by=12),seq(180,600,by=60))) { 
  
  upper.age <- ages
  mid.age <- 0.5*(c(0,ages[-length(ages)])+ages)
  avg.age <- var.age <- rep(NA,ncol(mat))	
  
  for (t in 1:ncol(mat)) { 
    pt <- cumsum(mat[,t]/sum(mat[,t]))
    opt <- optimize(f=findRate,interval=c(1e-9,10),ages=ages,prop=pt, maximum=TRUE)
     #plot(upper.age/12,pt, type="l")
     #points(upper.age/12,1-exp(-opt$maximum*ages), type="l",col=2)

    #print("show values")
    #print(findRate(par=0.1,ages=ages,prop=pt))
    #print(findRate(par=0.01,ages=ages,prop=pt))
    #print(findRate(par=opt$maximum,ages=ages,prop=pt))
    
    #print(opt)
    #get the median
    avg.age[t] <- 1/opt$maximum
    var.age[t] <- 1/(opt$maximum^2)		
  }
  
  return(list(avg.age=avg.age,var.age=var.age))
  
}



plotAge <- function(tmp,n.per.year=26,max.age=120) { 
  
  mat <- tmp$nt.store[tmp$i.inds,]
  mid.age <- 0.5*(c(0,tmp$age.classes[-length(tmp$age.classes)])+tmp$age.classes)
  avg.age <- var.age <- rep(NA,ncol(mat))	

  t.plot <- seq(1,ncol(mat),by=n.per.year)
  pt <- (t(mat)/colSums(mat))
  persp(t.plot[t.plot>tmp$burnin]/n.per.year,tmp$age.classes[tmp$age.classes<max.age]/12,
        pt[t.plot[t.plot>tmp$burnin],tmp$age.classes<max.age], 
        xlab="Time", ylab="Ages", zlab="Cases",  theta = 90, phi = 0)
  
  image(t.plot[t.plot>tmp$burnin]/n.per.year,tmp$age.classes[tmp$age.classes<max.age]/12,
        pt[t.plot[t.plot>tmp$burnin],tmp$age.classes<max.age], 
        xlab="Time", ylab="Ages")
  abline(v=tmp$vacc.t.multip/n.per.year)
  
  return(list(t.plot=t.plot[t.plot>tmp$burnin]/n.per.year,ages=tmp$age.classes[tmp$age.classes<max.age]/12, 
              pt=pt[t.plot[t.plot>tmp$burnin],tmp$age.classes<max.age]))
  
  }


tweakWAIFWtoMatchIncidence <- function(target.incidence=54, 
                                    test.multip=seq(25,40,length=20),
                                    age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)),  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                                    mort=c(rep(0,72),1),					#neglgible death till age 50 (could make more reaslitic)
                                    waifw=matrix(0.0002,73,73),					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                                    alpha=0.2,							#parameter governing sine wave seasonality in transmission	
                                    discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                                    waning.maternal.par=0.5,						#rate of waning of maternal immunity
                                    waning.immunity=0, 						#rate of waning of immunity in general
                                    vaccination=c(rep(0,68),rep(0.5,5)),				#vaccination strategy (e.g., this default: vaccinate at age 20, vaccinate 50% of the population)		
                                    time.step=1/4,							#time-step (1/4 of a month = 1 week)
                                    fert= c(rep(0,66),rep(0.0008,7)),					#pattern of fertility over age
                                    alpha.births = 0,						#pattern of seasonality in births	
                                    tmatch=140,burnin=10,Tmax=150, 
                                    do.plot=FALSE,	
                                    start.pop=c(),							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                                    start.pop.size=100000,						#starting total population size
                                    start.full.pop.struct=c(),					#starting full pop structure - for if you want to switch off vacc
                                    vacc.t.multip=52*10,						#time to multiply out the vacc by chosen value
                                    value.multip.vacc=1,						#value for 
                                    gamma=0.97){
  
  
  dist.target <- rep(NA,length(test.multip))
  
  for (k in 1:length(test.multip)){
  
  tmp <- iterateSIR(age.classes=age.classes,  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                    mort=mort,					#neglgible death till age 50 (could make more reaslitic)
                    waifw=waifw*test.multip[k],					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                    alpha=alpha,							#parameter governing sine wave seasonality in transmission	
                    discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                    waning.maternal.par=waning.maternal.par,						#rate of waning of maternal immunity
                    waning.immunity=waning.immunity, 						#rate of waning of immunity in general
                    vaccination=vaccination,				#vaccination strategy (e.g., this default: vaccinate at age 20, vaccinate 50% of the population)		
                    time.step=time.step,							#time-step (1/4 of a month = 1 week)
                    fert= fert,					#pattern of fertility over age
                    alpha.births = alpha.births,						#pattern of seasonality in births	
                    burnin=burnin,Tmax=Tmax, 
                    do.plot=FALSE,	
                    start.pop=start.pop,							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                    start.pop.size=start.pop.size,						#starting total population size
                    start.full.pop.struct=start.full.pop.struct,					#starting full pop structure - for if you want to switch off vacc
                    vacc.t.multip=vacc.t.multip,						#time to multiply out the vacc by chosen value
                    value.multip.vacc=value.multip.vacc,						#value for 
                    gamma=gamma)

  i.obtained <- sum(tmp$nt.store[tmp$i.inds,tmatch])
  dist.target[k] <- sum((i.obtained-target.incidence)^2)
  print(c(test.multip[k],i.obtained))
  
  }
  
    plot(test.multip, sqrt(dist.target), type="b",pch=19, xlab="Multiplier", ylab="Absolute distance", ylim=range(c(0,sqrt(dist.target))))
      best <- median(test.multip[which(dist.target==min(dist.target))])
      abline(v=best)

  return(list(test.multip=test.multip,dist.target=dist.target,waifw=waifw*best,i.obtained=i.obtained, best=best))
   
}



## Binary Search - this will only work if ~ linear 
tweakWAIFWtoMatchIncidenceBinary <- function(target.incidence=54, 
                                       test.multip=c(0.1,10),
                                       age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)),  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                                       mort=c(rep(0,72),1),					#neglgible death till age 50 (could make more reaslitic)
                                       waifw=matrix(0.0002,73,73),					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                                       alpha=0.2,							#parameter governing sine wave seasonality in transmission	
                                       discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                                       waning.maternal.par=0.5,						#rate of waning of maternal immunity
                                       waning.immunity=0, 						#rate of waning of immunity in general
                                       vaccination=c(rep(0,68),rep(0.5,5)),				#vaccination strategy (e.g., this default: vaccinate at age 20, vaccinate 50% of the population)		
                                       time.step=1/4,							#time-step (1/4 of a month = 1 week)
                                       fert= c(rep(0,66),rep(0.0008,7)),					#pattern of fertility over age
                                       alpha.births = 0,						#pattern of seasonality in births	
                                       tmatch=140,burnin=10,Tmax=150, 
                                       do.plot=FALSE,	
                                       start.pop=c(),							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                                       start.pop.size=100000,						#starting total population size
                                       start.full.pop.struct=c(),					#starting full pop structure - for if you want to switch off vacc
                                       vacc.t.multip=52*10,						#time to multiply out the vacc by chosen value
                                       value.multip.vacc=1,						#value for 
                                       gamma=0.97, step.size=0.1,tol=1e-5,maxtest=100){
  count <- 0; dist.from.target <- 100
  store.resp <- store.test <- rep(NA,maxtest) 
  while(count<maxtest | dist.from.target>tol){
    count <- count + 1
    new.test <- mean(test.multip)
  
    tmp <- iterateSIR(age.classes=age.classes,  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                      mort=mort,					   #neglgible death till age 50 (could make more reaslitic)
                      waifw=waifw*new.test,					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                      alpha=alpha,							#parameter governing sine wave seasonality in transmission	
                      discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                      waning.maternal.par=waning.maternal.par,						#rate of waning of maternal immunity
                      waning.immunity=waning.immunity, 						#rate of waning of immunity in general
                      vaccination=vaccination,				#vaccination strategy (e.g., this default: vaccinate at age 20, vaccinate 50% of the population)		
                      time.step=time.step,							#time-step (1/4 of a month = 1 week)
                      fert= fert,					#pattern of fertility over age
                      alpha.births = alpha.births,						#pattern of seasonality in births	
                      burnin=burnin,Tmax=Tmax, 
                      do.plot=FALSE,	
                      start.pop=start.pop,							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                      start.pop.size=start.pop.size,						#starting total population size
                      start.full.pop.struct=start.full.pop.struct,					#starting full pop structure - for if you want to switch off vacc
                      vacc.t.multip=vacc.t.multip,						#time to multiply out the vacc by chosen value
                      value.multip.vacc=value.multip.vacc,						#value for 
                      gamma=gamma)
    
    i.obtained <- sum(tmp$nt.store[tmp$i.inds,tmatch])
    
    dist.from.target <- abs(i.obtained-target.incidence)
    if (count==maxtest) break()
    
    if (i.obtained<target.incidence) test.multip[1] <- test.multip[1]+step.size
    if (i.obtained>target.incidence) test.multip[2] <- test.multip[2]-step.size
    store.resp[count] <- i.obtained
    store.test[count] <- new.test

    if (count>2) if(store.test[count]==store.test[count-2]) break()
    
    print(c(count,new.test,dist.from.target,test.multip,i.obtained))
  }
  
  store.test <- c(store.test,new.test)
  store.resp <- c(store.resp,i.obtained)
  
  plot(store.test, store.resp, type="b",pch=19, xlab="Multiplier", ylab="Incidence", ylim=range(c(store.resp, target.incidence),na.rm=TRUE))
  abline(h=target.incidence)
  abline(v=new.test)
  

  return(list(waifw=waifw*new.test,i.obtained=i.obtained, best=new.test,target.incidence=target.incidence))
  
}




tweakWAIFWtoMatchAge <- function(target.age=54, 
                                       test.multip=seq(25,40,length=20),
                                       age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)),  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                                       mort=c(rep(0,72),1),					#neglgible death till age 50 (could make more reaslitic)
                                       waifw=matrix(0.0002,73,73),					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                                       alpha=0.2,							#parameter governing sine wave seasonality in transmission	
                                       discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                                       waning.maternal.par=0.5,						#rate of waning of maternal immunity
                                       waning.immunity=0, 						#rate of waning of immunity in general
                                       vaccination=c(rep(0,68),rep(0.5,5)),				#vaccination strategy (e.g., this default: vaccinate at age 20, vaccinate 50% of the population)		
                                       time.step=1/4,							#time-step (1/4 of a month = 1 week)
                                       fert= c(rep(0,66),rep(0.0008,7)),					#pattern of fertility over age
                                       alpha.births = 0,						#pattern of seasonality in births	
                                       tmatch=140,burnin=10,Tmax=150, 
                                       do.plot=FALSE,	
                                       start.pop=c(),							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                                       start.pop.size=100000,						#starting total population size
                                       start.full.pop.struct=c(),					#starting full pop structure - for if you want to switch off vacc
                                       vacc.t.multip=52*10,						#time to multiply out the vacc by chosen value
                                       value.multip.vacc=1,						#value for 
                                       gamma=0.97){
  
  
  dist.target <- rep(NA,length(test.multip))
  
  for (k in 1:length(test.multip)){
    print(k)
    tmp <- iterateSIR(age.classes=age.classes,  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                      mort=mort,					#neglgible death till age 50 (could make more reaslitic)
                      waifw=waifw*test.multip[k],					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                      alpha=alpha,							#parameter governing sine wave seasonality in transmission	
                      discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                      waning.maternal.par=waning.maternal.par,						#rate of waning of maternal immunity
                      waning.immunity=waning.immunity, 						#rate of waning of immunity in general
                      vaccination=vaccination,				#vaccination strategy (e.g., this default: vaccinate at age 20, vaccinate 50% of the population)		
                      time.step=time.step,							#time-step (1/4 of a month = 1 week)
                      fert= fert,					#pattern of fertility over age
                      alpha.births = alpha.births,						#pattern of seasonality in births	
                      burnin=burnin,Tmax=Tmax, 
                      do.plot=FALSE,	
                      start.pop=start.pop,							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                      start.pop.size=start.pop.size,						#starting total population size
                      start.full.pop.struct=start.full.pop.struct,					#starting full pop structure - for if you want to switch off vacc
                      vacc.t.multip=vacc.t.multip,						#time to multiply out the vacc by chosen value
                      value.multip.vacc=value.multip.vacc,						#value for 
                      gamma=gamma)
    
    age.vals <- extractAge(tmp$nt.store[tmp$i.inds,],tmp$age.classes)
    
    age.obtained <- mean(age.vals$avg.age[tmatch],na.rm=TRUE)
    dist.target[k] <- sum((age.obtained-target.age)^2)
    print(age.obtained)
    
  }
  
  plot(test.multip, sqrt(dist.target), type="b",pch=19, xlab="Multiplier", ylab="Absolute distance", ylim=range(c(0,sqrt(dist.target)),na.rm=TRUE))
  best <- test.multip[which(dist.target==min(dist.target,na.rm=T))]
  abline(v=best)
  
  return(list(test.multip=test.multip,dist.target=dist.target,waifw=waifw*best,age.obtained=age.obtained,best=best))
  
}


tweakWAIFWtoMatchAgeBinary <- function(target.age=50, 
                                       test.multip=c(0.1,10),
                                       age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)),  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                                       mort=c(rep(0,72),1),					#neglgible death till age 50 (could make more reaslitic)
                                       waifw=matrix(0.0002,73,73),					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                                       alpha=0.2,							#parameter governing sine wave seasonality in transmission	
                                       discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                                       waning.maternal.par=0.5,						#rate of waning of maternal immunity
                                       waning.immunity=0, 						#rate of waning of immunity in general
                                       vaccination=c(rep(0,68),rep(0.5,5)),				#vaccination strategy (e.g., this default: vaccinate at age 20, vaccinate 50% of the population)		
                                       time.step=1/4,							#time-step (1/4 of a month = 1 week)
                                       fert= c(rep(0,66),rep(0.0008,7)),					#pattern of fertility over age
                                       alpha.births = 0,						#pattern of seasonality in births	
                                       tmatch=140,burnin=10,Tmax=150, 
                                       do.plot=FALSE,	
                                       start.pop=c(),							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                                       start.pop.size=100000,						#starting total population size
                                       start.full.pop.struct=c(),					#starting full pop structure - for if you want to switch off vacc
                                       vacc.t.multip=52*10,						#time to multiply out the vacc by chosen value
                                       value.multip.vacc=1,						#value for 
                                       gamma=0.97, step.size=0.1,tol=1e-3,maxtest=100){
  count <- 0; dist.from.target <- 100
  store.resp <- store.test <- rep(NA,maxtest) 
  while(count<maxtest | dist.from.target>tol){
    count <- count + 1
    new.test <- mean(test.multip)
    
    tmp <- iterateSIR(age.classes=age.classes,  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                      mort=mort,					   #neglgible death till age 50 (could make more reaslitic)
                      waifw=waifw*new.test,					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                      alpha=alpha,							#parameter governing sine wave seasonality in transmission	
                      discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                      waning.maternal.par=waning.maternal.par,						#rate of waning of maternal immunity
                      waning.immunity=waning.immunity, 						#rate of waning of immunity in general
                      vaccination=vaccination,				#vaccination strategy (e.g., this default: vaccinate at age 20, vaccinate 50% of the population)		
                      time.step=time.step,							#time-step (1/4 of a month = 1 week)
                      fert= fert,					#pattern of fertility over age
                      alpha.births = alpha.births,						#pattern of seasonality in births	
                      burnin=burnin,Tmax=Tmax, 
                      do.plot=FALSE,	
                      start.pop=start.pop,							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                      start.pop.size=start.pop.size,						#starting total population size
                      start.full.pop.struct=start.full.pop.struct,					#starting full pop structure - for if you want to switch off vacc
                      vacc.t.multip=vacc.t.multip,						#time to multiply out the vacc by chosen value
                      value.multip.vacc=value.multip.vacc,						#value for 
                      gamma=gamma)
    
    age.vals <- extractAge(tmp$nt.store[tmp$i.inds,],tmp$age.classes)
    age.obtained <- mean(age.vals$avg.age[tmatch],na.rm=TRUE)
    
    ## has to be an increase so do one over 
    dist.from.target <- abs((1/age.obtained)-(1/target.age))
    if (count==maxtest) break()
    
    if (is.na(age.obtained))
      
      if ((1/age.obtained)<(1/target.age)) test.multip[1] <- test.multip[1]+step.size
    if ((1/age.obtained)>(1/target.age)) test.multip[2] <- test.multip[2]-step.size
    store.resp[count] <- age.obtained
    store.test[count] <- new.test
    
    if (count>2) if(store.test[count]==store.test[count-2]) break()
    
    print(c(count,new.test,dist.from.target,test.multip))
  }
  
  store.test <- c(store.test,new.test)
  store.resp <- c(store.resp,age.obtained)
  
  plot(store.test, store.resp, type="b",pch=19, xlab="Multiplier", ylab="Age", ylim=range(c(store.resp, target.age),na.rm=TRUE))
  abline(h=target.age)
  abline(v=new.test)
  
  return(list(waifw=waifw*new.test,age.obtained=age.obtained, best=new.test,target.age=target.age))
  
}


tweakWAIFWtoMatchBoth <- function(target.incidence=54, target.age=50,
                                       test.multip=seq(25,40,length=20),
                                       age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)),  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                                       mort=c(rep(0,72),1),					#neglgible death till age 50 (could make more reaslitic)
                                       waifw=matrix(0.0002,73,73),				#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                                       alpha=0.2,							#parameter governing sine wave seasonality in transmission	
                                       discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                                       waning.maternal.par=0.5,						#rate of waning of maternal immunity
                                       waning.immunity=0, 						#rate of waning of immunity in general
                                       vaccination=c(rep(0,68),rep(0.5,5)),				#vaccination strategy (e.g., this default: vaccinate at age 20, vaccinate 50% of the population)		
                                       time.step=1/4,							#time-step (1/4 of a month = 1 week)
                                       fert= c(rep(0,66),rep(0.0008,7)),					#pattern of fertility over age
                                       alpha.births = 0,						#pattern of seasonality in births	
                                       tmatch=140,burnin=10,Tmax=150, 
                                       do.plot=FALSE,	
                                       start.pop=c(),							#starting population structur, # - could impose, otherwise takes the equ structure with imposed start # 
                                       start.pop.size=100000,						#starting total population size
                                       start.full.pop.struct=c(),					#starting full pop structure - for if you want to switch off vacc
                                       vacc.t.multip=52*10,						#time to multiply out the vacc by chosen value
                                       value.multip.vacc=1,
				       value.add.vacc=0,	#value for 
                                       gamma=0.97){
  
  
  i.obtained <- age.obtained <- dist.target <- dist.target.age <- rep(NA,length(test.multip))
  

  for (k in 1:length(test.multip)){
    
    tmp <- iterateSIR(age.classes=age.classes,  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                      mort=mort,					#neglgible death till age 50 (could make more reaslitic)
                      waifw=waifw*test.multip[k],					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                      alpha=alpha,							#parameter governing sine wave seasonality in transmission	
                      discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                      waning.maternal.par=waning.maternal.par,						#rate of waning of maternal immunity
                      waning.immunity=waning.immunity, 						#rate of waning of immunity in general
                      vaccination=vaccination,				#vaccination strategy (e.g., this default: vaccinate at age 20, vaccinate 50% of the population)		
                      time.step=time.step,							#time-step (1/4 of a month = 1 week)
                      fert= fert,					#pattern of fertility over age
                      alpha.births = alpha.births,						#pattern of seasonality in births	
                      burnin=burnin,Tmax=Tmax, 
                      do.plot=FALSE,	
                      start.pop=start.pop,							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                      start.pop.size=start.pop.size,						#starting total population size
                      start.full.pop.struct=start.full.pop.struct,					#starting full pop structure - for if you want to switch off vacc
                      vacc.t.multip=vacc.t.multip,						#time to multiply out the vacc by chosen value
                      value.multip.vacc=value.multip.vacc,						#value for 
                      value.add.vacc=value.add.vacc,gamma=gamma)
    
    i.obtained[k] <- sum(tmp$nt.store[tmp$i.inds,tmatch])
    dist.target[k] <- (i.obtained[k]-target.incidence)
    
    age.vals <- extractAge(tmp$nt.store[tmp$i.inds,],tmp$age.classes)
    age.obtained[k] <- mean(age.vals$avg.age[tmatch],na.rm=TRUE)
    dist.target.age[k] <- (age.obtained[k]-target.age)
    
    #print(c(test.multip[k],i.obtained[k], age.obtained[k]))
    
  }
  
  par(mfrow=c(2,1))
  plot(test.multip, (dist.target), type="b",pch=19, xlab="Multiplier", ylab="Absolute distance", ylim=range(c(0,(dist.target)),na.rm=TRUE))
  best <- median(test.multip[which(abs(dist.target)==min(abs(dist.target),na.rm=TRUE))])
  abline(v=best)
  abline(h=0)

  plot(test.multip, (dist.target.age), type="b",pch=19, xlab="Multiplier", ylab="Absolute distance age", ylim=range(c(0,(dist.target.age)),na.rm=TRUE))
  best.age <- median(test.multip[which(abs(dist.target.age)==min(abs(dist.target.age),na.rm=TRUE))])
  abline(v=best.age)
  abline(h=0)
  
  print("Multiplier")
  print(best)
  print(best.age)
  
  

  return(list(test.multip=test.multip,dist.target=dist.target,dist.target.age=dist.target.age,
              target.incidence=target.incidence,target.age=target.age,
              waifw=waifw*best,waifw.age=waifw*best.age,
              i.obtained=i.obtained,age.obtained=age.obtained, 
              best=best,best.age=best.age))
  
}




#function to get R0 from the waifw and Disease Free Eq
getR0 <- function(waifw,DFE){
  next.gen <- DFE*(1-exp(-waifw)); 
  next.gen <- pmin(next.gen,200000); #print(range(next.gen))
  cur.R0 <- Re(eigen(next.gen)$value[1])
  return(cur.R0)
}


extractStats <- function(incidence,avg.age,var.age, 
                            times.test=c(5,10,15)*52){
  ## time of peak   
  tpeak <- which(incidence==max(incidence))
  mats.rc <- cbind(incidence[times.test], avg.age[times.test], var.age[times.test])
  return((list(tpeak=tpeak,mats.rc=mats.rc)))
}


extractData <- function(){
  df <- read.csv("./data/QuarterlyEWMeaslesAge56to72.csv",header=TRUE, stringsAsFactors=FALSE)
  cases.male <- as.matrix(df[,seq(5,22,by=2)])
  cases.female <- as.matrix(df[,seq(6,22,by=2)])
  tot.cases <- cases.male+cases.female   
  upper.ages <- c(1,2,3,4,5,10,15,25)
  times <- df$Year+df$Quarter/4
  #image(times,upper.ages,log(tot.cases[,1:length(upper.ages)]), xlab="Year", ylab="Age")
  return(tot.cases=tot.cases)
}


sourceMeaslesUkAgeQuarterly <- function(){
  df <- read.csv("./data/QuarterlyEWMeaslesAge56to72.csv",header=TRUE, stringsAsFactors=FALSE)
  cases.male <- as.matrix(df[,seq(5,22,by=2)])
  cases.female <- as.matrix(df[,seq(6,22,by=2)])
  tot.cases <- cases.male+cases.female   
  upper.ages <- c(1,2,3,4,5,10,15,25)
  times <- df$Year+df$Quarter/4
  image(times,upper.ages,log(tot.cases[,1:length(upper.ages)]), xlab="Year", ylab="Age")
  abline(v=68, lwd=2)
  mid.ages <- 0.5*(c(0,upper.ages[-length(upper.ages)])+upper.ages)
  
  cols <- colorRampPalette(c("grey","red","blue"))(8)
  matplot(times,tot.cases[,1:length(upper.ages)], type="l", col=cols,lty=1, xlab="Year")
  abline(v=68, lwd=2)
  legend("topright", legend=c(upper.ages), col=cols, lty=1)

  prop <- tot.cases[,1:length(upper.ages)]/rowSums(tot.cases[,1:length(upper.ages)])
  cumprop <- t(apply(prop,1,cumsum))
  matplot(times,cumprop, type="l")
  
  avg.age <- rowSums(t(t(prop)*mid.ages))
  plot(times,avg.age, type="b", pch=19, xlab="", ylab="Average age")
  abline(v=68, lwd=2)
  
  mid.ages.adj <- mid.ages
  mid.ages.adj[6] <- 10 ## assume that all the cases in that wide age bin happen at the end to see how big it could be
  mid.ages.adj[7] <- 25 ## assume that all the cases in that wide age bin happen at the end to see how big it could be
  avg.age.adj <- rowSums(t(t(prop)*mid.ages.adj))
 
  plot(times,avg.age.adj, type="b", pch=19, xlab="", ylab="Average age", col=2, ylim=c(4,8))
  points(times,avg.age, type="b",col=1,pch=19)
  abline(v=68, lwd=2)
  ## assuming shift from the lower curve to the upper, looking at ~5-7 so around 2 years max  
  
  
  df1 <- read.csv("./data/QuarterlyEWMeaslesAge44to55.csv",header=TRUE)
  cases.male1 <- df[,seq(5,22,by=2)]
  cases.female1 <- df[,seq(6,22,by=2)]
  tot.cases1 <- cases.male1+cases.female1   
  upper.ages1 <- c(1,3,5,10,15,25)
  times1 <- df1$Year+df1$Quarter/4
  

}


getLikelihood <- function(par=c(log(45),log(0.4),log(129)),
                            idx.par=(1:3),start.par=c(45,0.4,129,5,0.2, 0.05,0.1),
                            age.classes=c(1:180,seq(240,600,by=60)),  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                            waifw=matrix(0.0002,187,187),
                            discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                            waning.maternal.par=0.5,						#rate of waning of maternal immunity
                            waning.immunity=0, 						#rate of waning of immunity in general
                            vaccination=c(rep(0,11),1,rep(0,175)),				#vaccination strategy (e.g., this default: vaccinate at 12 months, everyone)
                            time.step=1/2,							#time-step (1/4 of a month = 1 week)
                            mort = c(rep(0,186),1),
                            fert =  c(rep(0,186),1), ##note that if lambda is not 1; vast slippage in average age
                            alpha.births = 0,						#pattern of seasonality in births	
                            tmatch=140,burnin=10,Tmax=150, 
                            do.plot=FALSE,	
                            pop.n=45000000,start.pop=c(),							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                            start.pop.size=45000000,						#starting total population size
                            start.full.pop.struct=c(),					#starting full pop structure - for if you want to switch off vacc
                            vacc.t.multip=68*26,						#time to multiply out the vacc by chosen value
                            value.multip.vacc=1,						#value for 
                            gamma=0.97, data=extractData(), do.waifw=FALSE,show.result=FALSE){

  par.now <- start.par
  par.now[idx.par] <- exp(par)
  par.now[3] <- par.now[3]%%26 #13 if time.step = 1
  
  if (do.waifw) waifw <- get.smooth.WAIFW(age.class.boundries=age.classes/12, mu=par.now[4],sig=par.now[5], gam=par.now[6], delta=par.now[7])
  
 # print(range(waifw))
  
  multip <- par.now[1]
  alpha <- par.now[2]
  seas.offset <- par.now[3]
  
  #print(multip)
  #print(alpha)
  #print(seas.offset)
  
  n.per.yr <- 26
  nyrs <- 99
  nyrs.burnin <- 60
  
  ###Coverage gradually improved from approximately 50% during the 1970s to 86% when MMR vaccine replaced single antigen vaccine in 1988, and reached 92% in 1995.
  v.times <- rep(0,2600)
  j1 <- 68; j2<-88; j3 <- 92
  v.times[((j1-1)*26+1):((j2-1)*26+26)] <- seq(0.5,0.86,length=length(c(((j1-1)*26+1):((j2-1)*26+26))))
  v.times[((j2-1)*26+1):((j3-1)*26+26)] <- seq(0.86,0.92,length=length(c(((j2-1)*26+1):((j3-1)*26+26))))
  v.times[((j3-1)*26+26):2600] <- 0.92
    
  ## cut down on vaccination to avoid extinction (justifiable by efficacy < 100%)
  v.times <- v.times*0.9
  
  DFE <- c(rep(start.pop.size/(length(age.classes)-1),length(age.classes)-1),0)
  waifw <- scale.Waifw(par.now[1], DFE, waifw)
  waifw[!is.finite(waifw)] <- max(waifw[is.finite(waifw)]) ## tests suggests this is what is weird
  #print(range(DFE*(1-exp(-waifw))))
  #image(DFE*(1-exp(-waifw)))
    
  ## Initiate time-series
  tmp0 <- iterateSIR(age.classes=age.classes, 
                    mort=mort,	
                    waifw=waifw,
                    alpha=alpha,seas.offset=seas.offset,
                    waning.maternal.par=waning.maternal.par,waning.immunity=waning.immunity,
                    vaccination=vaccination,						
                    time.step=time.step,
                    fert= fert,
                    vacc.t.multip=vacc.t.multip,value.multip.vacc=value.multip.vacc,	#switch off vaccination from vacc.t.multip
                    burnin=(n.per.yr*nyrs.burnin),Tmax=(n.per.yr*nyrs)*2, 
                    do.plot=FALSE,
                    start.pop=DFE, #start where you scaled the WAIFW to? 
                    start.pop.size=pop.n,start.full.pop.struct=start.full.pop.struct,
                    v.times=rep(v.times*0,10),
                    gamma=gamma)
  
  
  #print(tmp0$nt.store[,tmp0$Tmax])
 # plot(colSums(tmp0$nt.store[tmp0$i.inds,(tmp0$Tmax-1000):tmp0$Tmax]), col=2, type="l")
 #  plot(colSums(tmp0$nt.store[,(tmp0$Tmax-1000):tmp0$Tmax]), col=2, type="l")
  print(c("ann.pop.growth",sum(tmp0$nt.store[,Tmax])/sum(tmp0$nt.store[,Tmax-26]))) #shold be ± 1.0015 for uK
  #print(sum(tmp0$nt.store[,tmp0$Tmax]))
  
  ## Run time-series
  tmp <- iterateSIR(age.classes=age.classes, 
                    mort=mort,	
                    waifw=waifw,
                    alpha=alpha,seas.offset=seas.offset,
                    waning.maternal.par=waning.maternal.par,waning.immunity=waning.immunity,
                    vaccination=vaccination,						
                    time.step=time.step,
                    fert= fert,
                    vacc.t.multip=vacc.t.multip,value.multip.vacc=value.multip.vacc,	#switch off vaccination from vacc.t.multip
                    burnin=0,#wwant no interactions in this part
                    Tmax=(n.per.yr*nyrs), 
                    do.plot=FALSE,
                    start.pop=c(),  
                    start.pop.size=c(),start.full.pop.struct=tmp0$nt.store[,tmp0$Tmax-26],#or -25? 
                    v.times=v.times,
                    gamma=gamma)
  
  # check connection 
  x <-26
  plot(c(colSums(tmp0$nt.store[tmp$i.inds,(tmp0$Tmax-50):(tmp0$Tmax-x)]),colSums(tmp$nt.store[tmp$i.inds,1:(26*3)])), col=2, type="l")
  abline(v=26*c(0:20)); abline(v=length((tmp0$Tmax-50):(tmp0$Tmax-x)),lwd=2)
  
  ## 2. Index 
  age1 <- tmp$age.classes<12
  age2 <- tmp$age.classes<24 & !age1
  age3 <- tmp$age.classes<36 & !age1 & !age2
  age4 <- tmp$age.classes<48 & !age1 & !age2 & !age3
  age5 <- tmp$age.classes<60 & !age1 & !age2 & !age3 & !age4
  age10 <- tmp$age.classes<120 & !age1 & !age2 & !age3 & !age4 & !age5 
  age15 <- tmp$age.classes<180 & !age1 & !age2 & !age3 & !age4 & !age5 & !age10
  age25 <- tmp$age.classes<300 & !age1 & !age2 & !age3 & !age4 & !age5 & !age10 & !age15
  age25more <- tmp$age.classes>=300
  
  tmatch <- rep(1:4,each=6)
  tmatch <- c(tmatch[1:6],2,tmatch[7:18],3,tmatch[19:24]) #pop the extras in the middle and at the end (check sensitivity to this)

  #print(tmp$i.inds[age15])
  #print(dim(tmp$nt.store))
  
  ## 3. Totals by age group over time 
  a1 <- colSums(tmp$nt.store[tmp$i.inds[age1],]) 
  a2 <- colSums(tmp$nt.store[tmp$i.inds[age2],]) 
  a3 <- colSums(tmp$nt.store[tmp$i.inds[age3],]) 
  a4 <- colSums(tmp$nt.store[tmp$i.inds[age4],]) 
  a5 <- colSums(tmp$nt.store[tmp$i.inds[age5],]) 
  a10 <- colSums(tmp$nt.store[tmp$i.inds[age10],]) 
  a15 <- colSums(tmp$nt.store[tmp$i.inds[age15],]) 
  #a15 <- tmp$nt.store[tmp$i.inds[age15],] ## this is just one dimensional under current age structure
  a25 <- colSums(tmp$nt.store[tmp$i.inds[age25],]) 
  a25more <- colSums(tmp$nt.store[tmp$i.inds[age25more],]) 
  
  #print("at start")
  #print(sum(tmp$nt.store[,56*26]))
  
  ## 4. Divide into quarters 
    
  rc <- matrix(NA,9,100*4); plt.yr <- c()
  for (j in 1:ceiling(tmp$Tmax/26)){
    this.year <- ((j-1)*26+1):((j-1)*26+26); 
    this.year.quarters <- ((j-1)*4+1):((j-1)*4+4); 
    
  # print(this.year)
  #  print(this.year.quarters)
    
    
    rc[1,this.year.quarters] <- c(sum(a1[this.year][tmatch==1]),sum(a1[this.year][tmatch==2]),
                                  sum(a1[this.year][tmatch==3]),sum(a1[this.year][tmatch==4]))
    rc[2,this.year.quarters] <- c(sum(a2[this.year][tmatch==1]),sum(a2[this.year][tmatch==2]),
                                  sum(a2[this.year][tmatch==3]),sum(a2[this.year][tmatch==4]))
    rc[3,this.year.quarters] <- c(sum(a3[this.year][tmatch==1]),sum(a3[this.year][tmatch==2]),
                                  sum(a3[this.year][tmatch==3]),sum(a3[this.year][tmatch==4]))
    rc[4,this.year.quarters] <- c(sum(a4[this.year][tmatch==1]),sum(a4[this.year][tmatch==2]),
                                  sum(a4[this.year][tmatch==3]),sum(a4[this.year][tmatch==4]))
    rc[5,this.year.quarters] <- c(sum(a5[this.year][tmatch==1]),sum(a5[this.year][tmatch==2]),
                                  sum(a5[this.year][tmatch==3]),sum(a5[this.year][tmatch==4]))
    rc[6,this.year.quarters] <- c(sum(a10[this.year][tmatch==1]),sum(a10[this.year][tmatch==2]),
                                  sum(a10[this.year][tmatch==3]),sum(a10[this.year][tmatch==4]))
    rc[7,this.year.quarters] <- c(sum(a15[this.year][tmatch==1]),sum(a15[this.year][tmatch==2]),
                                  sum(a15[this.year][tmatch==3]),sum(a15[this.year][tmatch==4]))
    rc[8,this.year.quarters] <- c(sum(a25[this.year][tmatch==1]),sum(a25[this.year][tmatch==2]),
                                  sum(a25[this.year][tmatch==3]),sum(a25[this.year][tmatch==4]))
    rc[9,this.year.quarters] <- c(sum(a25more[this.year][tmatch==1]),sum(a25more[this.year][tmatch==2]),
                                  sum(a25more[this.year][tmatch==3]),sum(a25more[this.year][tmatch==4]))
    plt.yr <- c(plt.yr,this.year.quarters)
  }                         
  
  ## pick out years from 56.25 to 76, which is the years of data
  plt.yr <- plt.yr[((56-1)*4+1):ncol(rc)][1:80]; #print(plt.yr)
  rc <- rc[,((56-1)*4+1):ncol(rc)]
  rc <- rc[,1:80]
  rc[is.na(rc)] <- 0
  
  #print(t(data[,1:8])-rc)
  
  #par(mfrow=c(2,2))
  matplot(plt.yr/4,t(rc), type="l", xlab="Time", ylab="Cases", col=c(1:5,rep("grey",4)),log="y", #put wide add range in grey
          lty=c(1,1,1,1,1,2,3,1,1), ylim=range(c(data+1,t(rc)+1),na.rm=TRUE)); # log="y",
  abline(v=50:80, lty=1,col="grey")
  matplot(plt.yr/4,data, add=TRUE, type="p",pch=c(rep(19,5),c(15,3,1,19)), col=c(1:5,rep("grey",4)))
  abline(v=50:80, lty=1,col="grey")
  
  
    #tst <- tmp$Tmax
  #plot(tmp$age.classes/12,tmp$nt.store[tmp$m.inds,tst]+tmp$nt.store[tmp$s.inds,tst]+
  #      tmp$nt.store[tmp$i.inds,tst]+tmp$nt.store[tmp$r.inds,tst], type="l")
  #plot(colSums(tmp$nt.store[tmp$s.inds,260:tst]), type="l", ylab="susc")
  #points(colSums(tmp0$nt.store[tmp0$s.inds,(tmp0$Tmax-tst):tmp0$Tmax]), col=2, type="l")
  #plot(colSums(tmp$nt.store[tmp$i.inds,26:tst]), type="l", ylab="infect")
  #points(colSums(tmp0$nt.store[tmp0$i.inds,(tmp0$Tmax-tst):tmp0$Tmax]), col=2, type="l")
  
  #print(range(colSums(tmp0$nt.store[tmp0$i.inds,(tmp0$Tmax-260):tmp0$Tmax])))
  
#  u <- (sum((t(data[,1:8])-rc)^2,na.rm=TRUE))
#  u <- log(sum((t(log(data[1:79,1:9]+1))-log(rc[,1:79]+1))^2,na.rm=TRUE))
#  u <- -sum(pmax(dpois(c(data[1:79,1:9]),c(t(rc[,1:79])),log=TRUE),-exp(300)))

  u <- -sum(pmax(dnbinom(c(data[1:79,1:9]),mu=c(t(rc[,1:79])),size=c(t(rc[,1:79])),log=TRUE),-exp(300)))
  
  #this one doubles the weight of that single age class one
#  u <- -sum(pmax(dpois(c(data[1:79,1:9]),c(t(rc[,1:79])),log=TRUE)*rep(c(rep(1,6),2,1,1),79),-exp(200)))
  
    print(u); print(exp(par))
  if (sum(rc)==0) {print("all zeros"); u <- exp(200)}
  
  #hist(log(c(data[1:79,1:8])/c(t(rc[,1:79]))))
  #print(mean(c(data[1:79,1:8])/c(t(rc[,1:79]))))
 
  if (!show.result) return(u) else return(list(rc=rc[,1:79],data=t(data[1:79,1:9]), tmp=tmp,u=u))

 }






getSyntheticLikelihood <- function(par=c(log(45),log(0.4),log(129)),
                          idx.par=(1:3),start.par=c(45,0.4,129,5,0.2, 0.05,0.1),
                          age.classes=c(1:180,seq(240,600,by=60)),  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                          waifw=matrix(0.0002,187,187),
                          discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                          waning.maternal.par=0.5,						#rate of waning of maternal immunity
                          waning.immunity=0, 						#rate of waning of immunity in general
                          vaccination=c(rep(0,11),1,rep(0,175)),				#vaccination strategy (e.g., this default: vaccinate at 12 months, everyone)
                          time.step=1/2,							#time-step (1/4 of a month = 1 week)
                          mort = c(rep(0,186),1),
                          fert =  c(rep(0,186),1), ##note that if lambda is not 1; vast slippage in average age
                          alpha.births = 0,						#pattern of seasonality in births	
                          tmatch=140,burnin=10,Tmax=150, 
                          do.plot=FALSE,	
                          pop.n=45000000,start.pop=c(),							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                          start.pop.size=45000000,						#starting total population size
                          start.full.pop.struct=c(),					#starting full pop structure - for if you want to switch off vacc
                          vacc.t.multip=68*26,						#time to multiply out the vacc by chosen value
                          value.multip.vacc=1,						#value for 
                          gamma=0.97, data=extractData(), do.waifw=FALSE,show.result=FALSE){
  
  par.now <- start.par
  par.now[idx.par] <- exp(par)
  par.now[3] <- par.now[3]%%26 #13 if time.step = 1
  
  if (do.waifw) waifw <- get.smooth.WAIFW(age.class.boundries=age.classes/12, mu=par.now[4],sig=par.now[5], gam=par.now[6], delta=par.now[7])
  
  # print(range(waifw))
  
  multip <- par.now[1]
  alpha <- par.now[2]
  seas.offset <- par.now[3]
  
  #print(multip)
  #print(alpha)
  #print(seas.offset)
  
  n.per.yr <- 26
  nyrs <- 99
  nyrs.burnin <- 60
  
  ###Coverage gradually improved from approximately 50% during the 1970s to 86% when MMR vaccine replaced single antigen vaccine in 1988, and reached 92% in 1995.
  v.times <- rep(0,2600)
  j1 <- 68; j2<-88; j3 <- 92
  v.times[((j1-1)*26+1):((j2-1)*26+26)] <- seq(0.5,0.86,length=length(c(((j1-1)*26+1):((j2-1)*26+26))))
  v.times[((j2-1)*26+1):((j3-1)*26+26)] <- seq(0.86,0.92,length=length(c(((j2-1)*26+1):((j3-1)*26+26))))
  v.times[((j3-1)*26+26):2600] <- 0.92
  
  ## cut down on vaccination to avoid extinction (justifiable by efficacy < 100%)
  v.times <- v.times*0.9
  
  DFE <- c(rep(start.pop.size/(length(age.classes)-1),length(age.classes)-1),0)
  waifw <- scale.Waifw(par.now[1], DFE, waifw)
  waifw[!is.finite(waifw)] <- max(waifw[is.finite(waifw)]) ## tests suggests this is what is weird
  #print(range(DFE*(1-exp(-waifw))))
  #image(DFE*(1-exp(-waifw)))
  
  ## Initiate time-series
  tmp0 <- iterateSIR(age.classes=age.classes, 
                     mort=mort,	
                     waifw=waifw,
                     alpha=alpha,seas.offset=seas.offset,
                     waning.maternal.par=waning.maternal.par,waning.immunity=waning.immunity,
                     vaccination=vaccination,						
                     time.step=time.step,
                     fert= fert,
                     vacc.t.multip=vacc.t.multip,value.multip.vacc=value.multip.vacc,	#switch off vaccination from vacc.t.multip
                     burnin=(n.per.yr*nyrs.burnin),Tmax=(n.per.yr*nyrs)*2, 
                     do.plot=FALSE,
                     start.pop=DFE, #start where you scaled the WAIFW to? 
                     start.pop.size=pop.n,start.full.pop.struct=start.full.pop.struct,
                     v.times=rep(v.times*0,10),
                     gamma=gamma)
  
  
  #print(tmp0$nt.store[,tmp0$Tmax])
  # plot(colSums(tmp0$nt.store[tmp0$i.inds,(tmp0$Tmax-1000):tmp0$Tmax]), col=2, type="l")
  #  plot(colSums(tmp0$nt.store[,(tmp0$Tmax-1000):tmp0$Tmax]), col=2, type="l")
  print(c("ann.pop.growth",sum(tmp0$nt.store[,Tmax])/sum(tmp0$nt.store[,Tmax-26]))) #shold be ± 1.0015 for uK
  #print(sum(tmp0$nt.store[,tmp0$Tmax]))
  
  ## Run time-series
  tmp <- iterateSIR(age.classes=age.classes, 
                    mort=mort,	
                    waifw=waifw,
                    alpha=alpha,seas.offset=seas.offset,
                    waning.maternal.par=waning.maternal.par,waning.immunity=waning.immunity,
                    vaccination=vaccination,						
                    time.step=time.step,
                    fert= fert,
                    vacc.t.multip=vacc.t.multip,value.multip.vacc=value.multip.vacc,	#switch off vaccination from vacc.t.multip
                    burnin=0,#wwant no interactions in this part
                    Tmax=(n.per.yr*nyrs), 
                    do.plot=FALSE,
                    start.pop=c(),  
                    start.pop.size=c(),start.full.pop.struct=tmp0$nt.store[,tmp0$Tmax-26],#or -25? 
                    v.times=v.times,
                    gamma=gamma)
  
  # check connection 
  x <-26
  plot(c(colSums(tmp0$nt.store[tmp$i.inds,(tmp0$Tmax-50):(tmp0$Tmax-x)]),colSums(tmp$nt.store[tmp$i.inds,1:(26*3)])), col=2, type="l")
  abline(v=26*c(0:20)); abline(v=length((tmp0$Tmax-50):(tmp0$Tmax-x)),lwd=2)
  
  ## 2. Index 
  age1 <- tmp$age.classes<12
  age2 <- tmp$age.classes<24 & !age1
  age3 <- tmp$age.classes<36 & !age1 & !age2
  age4 <- tmp$age.classes<48 & !age1 & !age2 & !age3
  age5 <- tmp$age.classes<60 & !age1 & !age2 & !age3 & !age4
  age10 <- tmp$age.classes<120 & !age1 & !age2 & !age3 & !age4 & !age5 
  age15 <- tmp$age.classes<180 & !age1 & !age2 & !age3 & !age4 & !age5 & !age10
  age25 <- tmp$age.classes<300 & !age1 & !age2 & !age3 & !age4 & !age5 & !age10 & !age15
  age25more <- tmp$age.classes>=300
  
  tmatch <- rep(1:4,each=6)
  tmatch <- c(tmatch[1:6],2,tmatch[7:18],3,tmatch[19:24]) #pop the extras in the middle and at the end (check sensitivity to this)
  
  #print(tmp$i.inds[age15])
  #print(dim(tmp$nt.store))
  
  ## 3. Totals by age group over time 
  a1 <- colSums(tmp$nt.store[tmp$i.inds[age1],]) 
  a2 <- colSums(tmp$nt.store[tmp$i.inds[age2],]) 
  a3 <- colSums(tmp$nt.store[tmp$i.inds[age3],]) 
  a4 <- colSums(tmp$nt.store[tmp$i.inds[age4],]) 
  a5 <- colSums(tmp$nt.store[tmp$i.inds[age5],]) 
  a10 <- colSums(tmp$nt.store[tmp$i.inds[age10],]) 
  a15 <- colSums(tmp$nt.store[tmp$i.inds[age15],]) 
  #a15 <- tmp$nt.store[tmp$i.inds[age15],] ## this is just one dimensional under current age structure
  a25 <- colSums(tmp$nt.store[tmp$i.inds[age25],]) 
  a25more <- colSums(tmp$nt.store[tmp$i.inds[age25more],]) 
  
  #print("at start")
  #print(sum(tmp$nt.store[,56*26]))
  
  ## 4. Divide into quarters 
  
  rc <- matrix(NA,9,100*4); plt.yr <- c()
  for (j in 1:ceiling(tmp$Tmax/26)){
    this.year <- ((j-1)*26+1):((j-1)*26+26); 
    this.year.quarters <- ((j-1)*4+1):((j-1)*4+4); 
    
    # print(this.year)
    #  print(this.year.quarters)
    
    
    rc[1,this.year.quarters] <- c(sum(a1[this.year][tmatch==1]),sum(a1[this.year][tmatch==2]),
                                  sum(a1[this.year][tmatch==3]),sum(a1[this.year][tmatch==4]))
    rc[2,this.year.quarters] <- c(sum(a2[this.year][tmatch==1]),sum(a2[this.year][tmatch==2]),
                                  sum(a2[this.year][tmatch==3]),sum(a2[this.year][tmatch==4]))
    rc[3,this.year.quarters] <- c(sum(a3[this.year][tmatch==1]),sum(a3[this.year][tmatch==2]),
                                  sum(a3[this.year][tmatch==3]),sum(a3[this.year][tmatch==4]))
    rc[4,this.year.quarters] <- c(sum(a4[this.year][tmatch==1]),sum(a4[this.year][tmatch==2]),
                                  sum(a4[this.year][tmatch==3]),sum(a4[this.year][tmatch==4]))
    rc[5,this.year.quarters] <- c(sum(a5[this.year][tmatch==1]),sum(a5[this.year][tmatch==2]),
                                  sum(a5[this.year][tmatch==3]),sum(a5[this.year][tmatch==4]))
    rc[6,this.year.quarters] <- c(sum(a10[this.year][tmatch==1]),sum(a10[this.year][tmatch==2]),
                                  sum(a10[this.year][tmatch==3]),sum(a10[this.year][tmatch==4]))
    rc[7,this.year.quarters] <- c(sum(a15[this.year][tmatch==1]),sum(a15[this.year][tmatch==2]),
                                  sum(a15[this.year][tmatch==3]),sum(a15[this.year][tmatch==4]))
    rc[8,this.year.quarters] <- c(sum(a25[this.year][tmatch==1]),sum(a25[this.year][tmatch==2]),
                                  sum(a25[this.year][tmatch==3]),sum(a25[this.year][tmatch==4]))
    rc[9,this.year.quarters] <- c(sum(a25more[this.year][tmatch==1]),sum(a25more[this.year][tmatch==2]),
                                  sum(a25more[this.year][tmatch==3]),sum(a25more[this.year][tmatch==4]))
    plt.yr <- c(plt.yr,this.year.quarters)
  }                         
  
  ## pick out years from 56.25 to 76, which is the years of data
  plt.yr <- plt.yr[((56-1)*4+1):ncol(rc)][1:80]; #print(plt.yr)
  rc <- rc[,((56-1)*4+1):ncol(rc)]
  rc <- rc[,1:80]
  rc[is.na(rc)] <- 0
  
  #print(t(data[,1:8])-rc)
  
  #par(mfrow=c(2,2))
  matplot(plt.yr/4,t(rc), type="l", xlab="Time", ylab="Cases", col=c(1:5,rep("grey",4)),log="y", #put wide add range in grey
          lty=c(1,1,1,1,1,2,3,1,1), ylim=range(c(data+1,t(rc)+1),na.rm=TRUE)); # log="y",
  abline(v=50:80, lty=1,col="grey")
  matplot(plt.yr/4,data, add=TRUE, type="p",pch=c(rep(19,5),c(15,3,1,19)), col=c(1:5,rep("grey",4)))
  abline(v=50:80, lty=1,col="grey")
  
  
  #mean distance
  
  u <- sum((apply(rc,1,mean,na.rm=TRUE)-apply(data,2,mean,na.rm=TRUE))^2)+
    sum((apply(rc,1,quantile,0.95,na.rm=TRUE)-apply(data,2,quantile,0.95,na.rm=TRUE))^2)+
    sum((apply(rc,1,quantile,0.05,na.rm=TRUE)-apply(data,2,quantile,0.05,na.rm=TRUE))^2)
    

  sim1 <- ts(data = colSums(rc,na.rm=TRUE), start = c(1,1), deltat = 1/4)
  a1 <- spec.ar(sim1,plot=F);
  
  data1 <- ts(data = rowSums(data,na.rm=TRUE), start = c(1,1),  deltat = 1/4)
  a2 <- spec.ar(data1,plot=F);
  
#print(median(a1$freq[a1$spec==max(a1$spec)]))
#print(median(a2$freq[a2$spec==max(a2$spec)]))

  u <- u+10*((median(a1$freq[a1$spec==max(a1$spec)])-median(a2$freq[a2$spec==max(a2$spec)]))^2)

  print(u); print(exp(par))
  if (sum(rc)==0) {print("all zeros"); u <- exp(200)}
  
  if (!show.result) return(u) else return(list(rc=rc[,1:79],data=t(data[1:79,1:9]), tmp=tmp,u=u))
  
}



exploreUnity <- function(){

 require(tsiR)

  n.per.year <- 26
  time.step <- 1/2
  
  #starting parameters
  nyrs <- 70
  nyrs.burnin <- 35
  nyrs.match <- 45
  vacc.t.multip <- 50*n.per.year 
  pop.n <- 500000
  vacc.cover <- 0
  value.multip.vacc <- 0
  
  age.classes <- c(1,5,10,50)*12
  mort <- c(rep(0,3),1)
  fert <-  c(rep(0,3),rep(1,1))
  waning.maternal <- 0.5
  
  ## tweak this to adjust the rate of increase to get something close to equilbirium (corresponds to lambda=1)
  xx <- findStableStruct(age.classes=age.classes, 
                         mort=mort,fert =fert,time.step = time.step)
  xx$lambda
  
  scaled.flat.waifw <- scale.Waifw(R0=25, DFE=xx$stable.age*pop.n, waifw=matrix(0.0002,length(age.classes),length(age.classes)))
  
  ## 1. Flat Waifw
  tmp <- iterateSIR(age.classes=age.classes, mort=mort,	
                    waifw=scaled.flat.waifw,
                    alpha=0.2,
                    waning.maternal.par=waning.maternal,waning.immunity=0,
                    vaccination=c(0,0,0,0),						
                    time.step=time.step,
                    fert= fert,
                    vacc.t.multip=vacc.t.multip,value.multip.vacc=value.multip.vacc,	#switch off vaccination from vacc.t.multip
                    burnin=(n.per.year*nyrs.burnin),Tmax=(n.per.year*nyrs), 
                    do.plot=F,start.pop=c(),start.pop.size=pop.n,start.full.pop.struct=c(),gamma=0.97, do.unity=TRUE)
  
   age.vals <- extractAge(tmp$nt.store[tmp$i.inds,],tmp$age.classes)
 

   #make a datafrmame, subsmapling the cases with reporting rate 
   dataSim <- data.frame(time=1:ncol(tmp$nt.store), 
		susc=colSums(tmp$nt.store[tmp$s.inds,]),
		cases=rbinom(ncol(tmp$nt.store),floor(colSums(tmp$nt.store[tmp$i.inds,])),0.5),
		births=tmp$store.births, pop=colSums(tmp$nt.store), unity.estimate=tmp$store.beta.unity)
   dataSim <- dataSim[1000:nrow(dataSim),]

   dataSim$time <- dataSim$time-dataSim$time[1]+1

  sbar <- mean(dataSim$susc/dataSim$pop)
  runRes <- runtsir(data= dataSim, IP = 2, xreg = "cumcases", regtype="gaussian",
         	alpha = 0.97, sbar = sbar,family = "gaussian", pred="step-ahead",
		link = "identity",method = "negbin", nsim = 10)

  par(mfrow=c(1,3))
  plot(runRes$beta, type="b", ylim=range(c(runRes$beta,dataSim$unity.estimate)), xlab="", ylab=expression(beta))
  points(rep(c(26,1:25),100)[1:nrow(dataSim)],dataSim$unity.estimate, type="b", col="red")
  legend("topleft",legend=c("tsir","unity"), col=c("black","red"), bty="n",lty=1)

   subst <- 1:(26*4)   
   plot(rowSums(runRes$simS)[subst], ylim=range(c(rowSums(runRes$simS)[subst],dataSim$susc[subst])), xlab="time", ylab="susceptibles", type="b")
   points(dataSim$susc, type="b",col="red") 

   plot((rowSums(runRes$simS)[subst]-mean(rowSums(runRes$simS)))/sd(rowSums(runRes$simS)), xlab="time", ylab="susceptibles", type="b")
   points((dataSim$susc-mean(dataSim$susc))/sd(dataSim$susc), type="b",col="red") 


}



#### Do some messing around  ############################################################################################################
####

    # Coverage gradually improved from approximately 50% during the 1970s to 86% when MMR vaccine replaced single antigen vaccine in 1988, 
    # and reached 92% in 1995.
 

smallTest <- function(){

  n.per.year <- 26
  time.step <- 1/2
  
  #starting parameters
  nyrs <- 70
  nyrs.burnin <- 35
  nyrs.match <- 45
  vacc.t.multip <- 50*n.per.year 
  pop.n <- 500000
  vacc.cover <- 0
  value.multip.vacc <- 0
  
  
  age.classes <- c(1,5,10,50)*12
  mort <- c(rep(0,3),1)
  fert <-  c(rep(0,3),rep(1,1))
  waning.maternal <- 0.99
  
  ## tweak this to adjust the rate of increase to get something close to equilbirium (corresponds to lambda=1)
  xx <- findStableStruct(age.classes=age.classes, 
                         mort=mort,fert =fert,time.step = time.step)
  xx$lambda
  
  scaled.flat.waifw <- scale.Waifw(R0=R0, DFE=xx$stable.age*pop.n, waifw=matrix(0.0002,length(age.classes),length(age.classes)))
  
  ## 1. Flat Waifw
  tmp <- iterateSIR(age.classes=age.classes, mort=mort,	
                    waifw=scaled.flat.waifw,
                    alpha=0,
                    waning.maternal.par=waning.maternal,waning.immunity=0,
                    vaccination=c(0,0,0,0),						
                    time.step=time.step,
                    fert= fert,
                    vacc.t.multip=vacc.t.multip,value.multip.vacc=value.multip.vacc,	#switch off vaccination from vacc.t.multip
                    burnin=(n.per.year*nyrs.burnin),Tmax=(n.per.year*nyrs), 
                    do.plot=F,start.pop=c(),start.pop.size=pop.n,start.full.pop.struct=c(),gamma=0.97)
  
  age.vals <- extractAge(tmp$nt.store[tmp$i.inds,],tmp$age.classes)
  
  
}


exploreResults <- function(){

  time.step <- 0.5
  
	## tweak this to adjust the rate of increase to get something close to equilbirium (corresponds to lambda=1)
	xx <- findStableStruct(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
				mort= c(rep(0,72),1),fert =  c(rep(0,66),rep(0.0007,7)), time.step = time.step)
	plot(xx$age.classes/12,xx$stable.age,xlab="age", ylab="", type="l")
	xx$lambda
	## I have picked out the survival / fertility I got from this for the examples that follow 
	## easier to control/understand if pop isn't growing like crazy, or shrinking. 

	## check out the time-series - set very low vaccination 
	## with a flat waifw, putting beta into the waifw matrix at all points
	tmp <- iterateSIR(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), mort=c(rep(0,72),1),	
                          waifw=matrix(0.0002,73,73),alpha=0.2,
                         waning.maternal.par=0.5,waning.immunity=0,
			vaccination=c(rep(0,68),rep(0.001,5)),time.step=time.step,
                         fert= c(rep(0,66),rep(0.0007,7)),
			burnin=(52*20),Tmax=(52*30), do.plot=TRUE,start.pop=c(),gamma=0.97)

	## now with lower seasonality to show the annual (rather than biennial) dynamics
	tmp1 <- iterateSIR(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), mort=c(rep(0,72),1),	
                          waifw=matrix(0.0002,73,73),alpha=0.1,
                         waning.maternal.par=0.5,waning.immunity=0,
			vaccination=c(rep(0,68),rep(0.001,5)),time.step=time.step,
                         fert= c(rep(0,66),rep(0.0007,7)),
			burnin=(52*20),Tmax=(52*30), do.plot=TRUE,start.pop=c(),gamma=0.97)


	par(mfrow=c(1,2), bty="l")
  	plot((tmp$burnin:tmp$Tmax)/52, colSums(tmp$nt.store[tmp1$i.inds,tmp1$burnin:tmp$Tmax]), type="l",col=2,xlab="Time", ylab="Total cases")
	points((tmp$burnin:tmp$Tmax)/52, colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]), type="l",col=1)
	abline(v=c(0:100),lty=1,col="grey")

	plot(tmp1$age.classes,tmp1$nt.store[tmp$i.inds,52*30]/sum(tmp1$nt.store[tmp$i.inds,52*30]), type="l", xlim=c(0,60), xlab="Age (months)", ylab="Relative propn infected")
	points(tmp$age.classes,tmp$nt.store[tmp$i.inds,52*30]/sum(tmp$nt.store[tmp$i.inds,52*30]),type="l",col=2)
	legend("topright",legend=c("annual dynamics", "beinnial dynamics"),lty=1,col=c(1,2), bty="n")



	##look across an increasing magnitude of vaccination (shown in colours - red is no vaccination, blue is lots) 
	## this function is a wrapper to the previous - so letting fertility mort, etc be at their defaults. 
	## have also implemented some waning of immunity, because otherwise vaccination is entirely uninteresting in its effects, 
	## since all mothers are immunized. 
	a1<-exploreVaccination(waifw=matrix(0.002,73,73), 
				waning.maternal.par=0.2,
				waning.immunity=0.005,n.test=10,vaccination.template=c(rep(0,66),rep(1,7)), do.plot=TRUE)

	par(mfrow=c(1,3))
	plot(a1$vacc.test,a1$avg.age.infection, type="b",pch=19,xlab="Vaccination coverage", ylab="Average age of infection (years)") 
	plot(a1$vacc.test,a1$age.peak.infected, type="b",pch=19,xlab="Vaccination coverage", ylab="Age peak infection (years)")
	plot(a1$vacc.test,a1$number.infected, type="b",pch=19,xlab="Vaccination coverage", ylab="Number infected per year")
	


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
			alpha=0.2,
			discrete.seas=outall$sea_beta,
                         waning.maternal.par=0.5,waning.immunity=0,
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
                         waning.maternal.par=0.5,waning.immunity=0,
			vaccination=c(rep(0,68),rep(0.001,5)),						
			time.step=1/4,
                         fert=  c(rep(0,66),rep(0.1,7))*chk,
			burnin=(52*20),Tmax=(52*30), do.plot=TRUE,start.pop=c(),start.pop.size=outall$pop[1],gamma=0.97)


	par(mfrow=c(1,3))
	chs <- dfcases$state=="Florida"
	plot(times.cases[chs],(dfcases$percent_specimen_positive[chs])/max(dfcases$percent_specimen_positive[chs]),type="b", xlab="", ylab="cases", lty=1,lwd=2,pch=19)
 	points(((tmp$burnin:tmp$Tmax)-tmp$burnin)/52, colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax])/max(colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax])), type="l", col="blue")
	abline(v=c(0:100),lty=3)
	title("Florida")
		
	chs <- dfcases$state=="Texas"
	plot(times.cases[chs],(dfcases$percent_specimen_positive[chs])/max(dfcases$percent_specimen_positive[chs]),type="b", xlab="", ylab="cases", lty=1,lwd=2,pch=19)
	points(((tmp$burnin:tmp$Tmax)-tmp$burnin)/52, colSums(tmp1$nt.store[tmp1$i.inds,tmp$burnin:tmp$Tmax])/max(colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax])), 
		type="l", col="blue")
	abline(v=c(0:100),lty=3)
	title("Texas")
	
	#plot(tmp1$age.classes,tmp1$nt.store[tmp$i.inds,52*30]/sum(tmp1$nt.store[tmp$i.inds,52*30]), type="l", xlim=c(0,60),xlab="Age (months)", ylab="Relative propn infected", ylim=c(0,0.01))
	#points(tmp$age.classes,tmp$nt.store[tmp$i.inds,52*30]/sum(tmp$nt.store[tmp$i.inds,52*30]),type="l",col=2)

	indexes <- seq((52*29),(52*30),by=4)
	cols1 <- colorRampPalette(c(rgb(0,0,1,1), rgb(0,0,1,0)), alpha = TRUE)(7); cols1 <- c(cols1,cols1[6:1])
	cols2 <- colorRampPalette(c(rgb(1,0,0,1), rgb(1,0,0,0)), alpha = TRUE)(7); cols2 <- c(cols2,cols2[6:1])

	matplot(tmp1$age.classes,tmp1$nt.store[tmp$i.inds,indexes]/sum(tmp1$nt.store[tmp$i.inds,52*30]), type="l",lty=1,
		xlim=c(0,58), xlab="Age (months)", ylab="Relative propn infected", ylim=c(0,0.02),col=cols2)
	matplot(tmp$age.classes,tmp$nt.store[tmp$i.inds,indexes]/sum(tmp$nt.store[tmp$i.inds,52*30]),type="l",lty=1,col=cols1, add=TRUE)
	legend("topright",legend=c("Florida", "Texas"),lty=1,col=c(4,2), bty="n")

	points(tmp$age.classes,tmp1$nt.store[tmp$i.inds,(52-16)*30]/sum(tmp1$nt.store[tmp$i.inds,(52-16)*30]),type="l",col=1,lwd=2,lty=3)
	points(tmp$age.classes,tmp$nt.store[tmp$i.inds,(52-16)*30]/sum(tmp$nt.store[tmp$i.inds,(52-16)*30]),type="l",col=1,lwd=2,lty=3)


}




exploreFunkyWAIFW <- function(){

	fname <- ""
	load(paste(fname,"data/prelimRSVstatefitFlorida.RData",sep=""))

	time.step <- 1/2
	
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
				mort= c(rep(0,72),1),fert =  c(rep(0,66),rep(1,7)), time.step=time.step)

	# figure out what fertility should be to get these weekly births 
	chk <- findFert(target.births=mean(outall$births),fert.multiplier=seq(0.00001,0.002,length=10),popsize=outall$pop[1],
			age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
			mort=c(rep(1e-9,72),1), fert =  c(rep(0,66),rep(1,7)), time.step=time.step)

	#get stable age structure with chk in place
	xx <- findStableStruct(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), 
				mort= c(rep(0,72),1),fert =  c(rep(0,66),rep(1,7))*chk, time.step=time.step)

	## rescale the waifw
	scaled.waifw.mat <- scale.Waifw(R0=median(outall$sea_beta)*outall$pop[1], DFE=xx$stable.age*outall$pop[1], waifw=waifw.mat)

	tmp <- iterateSIR(age.classes=c(1:60,seq(72,120,by=12),seq(180,600,by=60)), mort=c(rep(0,72),1),	
                         waifw=scaled.waifw.mat,alpha=0.2,
			                    discrete.seas=outall$sea_beta,
                         waning.maternal.par=0.5,waning.immunity=0,
			                  vaccination=c(rep(0,68),rep(0.001,5)),						
			                    time.step=time.step,
                         fert=  c(rep(0,66),rep(1,7))*chk,
			                    burnin=(52*20),Tmax=(52*30), do.plot=TRUE,start.pop=c(),start.pop.size=outall$pop[1],gamma=0.97)


	dfcases <- read.csv(paste(fname,"data/RSV_data_US.csv",sep=""))
	times.cases <- dfcases$year+dfcases$week/52-2020+4

	par(mfrow=c(1,3))
	chs <- dfcases$state=="Florida"
	plot(times.cases[chs],(dfcases$percent_specimen_positive[chs])/max(dfcases$percent_specimen_positive[chs]),type="b", xlab="", ylab="cases", lty=1,lwd=2,pch=19)
 	points(((tmp$burnin:tmp$Tmax)-tmp$burnin)/52, colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax])/max(colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax])), 
		type="l", col="blue")
	abline(v=c(0:100),lty=3)
	title("Florida")



}




### wrapper model, also adjusted to be measles like parameters, aka, time-step =1/2, and 26 where there is 52

aFullComparison <- function(R0=30, remove.vacc= TRUE, 
                            waifw.test=get.smooth.WAIFW(age.class.boundries=c(1:120,seq(180,600,by=60))/12, mu=30,sig=1, gam=0.08, delta=0.01), #diagonal
                            #waifw.mat=get.smooth.WAIFW(age.class.boundries=c(1:120,seq(180,600,by=60))/12, mu=6,sig=0.4, gam=0.07, delta=0.001) ##blobby
                            full.plot=TRUE){
  
  n.per.year <- 26
  time.step <- 1/2
  
  #starting parameters
  nyrs <- 45#50
  nyrs.burnin <- 25
  nyrs.match <- 35
  vacc.t.multip <- 35*n.per.year 
  pop.n <- 100000
  ## remove vacc
  if (remove.vacc) { 
    vacc.cover <- 0.8
    value.multip.vacc <- 0 } else {
  ## introduce vacc
    vacc.cover <- 0.0008
    value.multip.vacc <- 1000
  }
  age.classes <- c(1:120,seq(180,600,by=60))
  mort <- c(rep(0,127),1)
  fert <-  c(rep(0,121),rep(0.000695,7))
  
  ## tweak this to adjust the rate of increase to get something close to equilbirium (corresponds to lambda=1)
  xx <- findStableStruct(age.classes=age.classes, 
                         mort=mort,fert =fert)
  #xx$lambda
  scaled.flat.waifw <- scale.Waifw(R0=R0, DFE=xx$stable.age*pop.n, waifw=matrix(0.0002,length(age.classes),length(age.classes)))
  
  ## 1. Flat Waifw
  tmp <- iterateSIR(age.classes=age.classes, mort=mort,	
                    waifw=scaled.flat.waifw,
                    alpha=0,
                    waning.maternal.par=0.5,waning.immunity=0,
                    vaccination=c(rep(0,11),vacc.cover,rep(0,116)),						
                    time.step=time.step,
                    fert= fert,
                    vacc.t.multip=vacc.t.multip,value.multip.vacc=value.multip.vacc,	#switch off vaccination from vacc.t.multip
                    burnin=(n.per.year*nyrs.burnin),Tmax=(n.per.year*nyrs), 
                    do.plot=F,start.pop=c(),start.pop.size=pop.n,start.full.pop.struct=c(),gamma=0.97)
  
  age.vals <- extractAge(tmp$nt.store[tmp$i.inds,],tmp$age.classes)
  
  
  ## 2. THEORETICAL - very focussed by age
  #waifw.mat <- get.smooth.WAIFW(age.class.boundries=age.classes/12, mu=6,sig=0.4, gam=0.07, delta=0.001)
  scaled.waifw.mat <- scale.Waifw(R0=R0, DFE=xx$stable.age*pop.n, waifw=waifw.test)

  ## match incidence to the previous 
  adj <- tweakWAIFWtoMatchIncidence(target.incidence=sum(tmp$nt.store[tmp$i.inds,(nyrs.match*n.per.year):((nyrs.match+1)*n.per.year)]),
                                    test.multip=seq(0.5,10,length=25),  
                                    age.classes=tmp$age.classes,  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                                    mort=mort,					#neglgible death till age 50 (could make more reaslitic)
                                    waifw=scaled.waifw.mat,					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                                    alpha=0,							#parameter governing sine wave seasonality in transmission	
                                    discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                                    waning.maternal.par=0.5,						#rate of waning of maternal immunity
                                    waning.immunity=0, 						#rate of waning of immunity in general
                                    vaccination=c(rep(0,11),vacc.cover,rep(0,116)),						
                                    time.step=time.step,							#time-step (1/4 of a month = 1 week)
                                    fert=fert,					#pattern of fertility over age
                                    alpha.births = 0,						#pattern of seasonality in births	
                                    tmatch=(nyrs.match*n.per.year):((nyrs.match+1)*n.per.year),burnin=500,Tmax=(nyrs.match+3)*n.per.year, 
                                    do.plot=FALSE,	
                                    start.pop=c(),							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                                    start.pop.size=pop.n,						#starting total population size
                                    start.full.pop.struct=c(),					#starting full pop structure - for if you want to switch off vacc
                                    vacc.t.multip=100*100,						#time to multiply out the vacc by chosen value - set to far beyond Tmax!!
                                    value.multip.vacc=1,						#value for 
                                    gamma=0.97)#, tol=0.1,max.test=30)
  
  if(adj$best==1) {new.range <- seq(0.01,1.1,length=20); ss <- 0.01 }
  if(adj$best==5) {new.range <- seq(4.5,20,length=20); ss <- 0.25 }
  if(adj$best!=5 & adj$best!=1) {new.range <- seq(adj$best-diff(adj$test.multip)[1],adj$best+diff(adj$test.multip)[1],length=20); ss <- 0.01}
  first.best.incidence <- adj$best
  
  #adj <- tweakWAIFWtoMatchIncidenceBinary(target.incidence=sum(tmp$nt.store[tmp$i.inds,(nyrs.match*n.per.year):((nyrs.match+1)*n.per.year)]),
  adj <- tweakWAIFWtoMatchIncidence(target.incidence=sum(tmp$nt.store[tmp$i.inds,(nyrs.match*n.per.year):((nyrs.match+1)*n.per.year)]),
                                    test.multip=new.range,#[c(1,20)],  
                                    age.classes=tmp$age.classes,  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                                    mort=mort,					#neglgible death till age 50 (could make more reaslitic)
                                    waifw=scaled.waifw.mat,					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                                    alpha=0,							#parameter governing sine wave seasonality in transmission	
                                    discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                                    waning.maternal.par=0.5,						#rate of waning of maternal immunity
                                    waning.immunity=0, 						#rate of waning of immunity in general
                                    vaccination=c(rep(0,11),vacc.cover,rep(0,116)),						
                                    time.step=time.step,							#time-step (1/4 of a month = 1 week)
                                    fert=fert,					#pattern of fertility over age
                                    alpha.births = 0,						#pattern of seasonality in births	
                                    tmatch=(nyrs.match*n.per.year):((nyrs.match+1)*n.per.year),burnin=500,Tmax=(nyrs.match+3)*n.per.year, 
                                    do.plot=FALSE,	
                                    start.pop=c(),							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                                    start.pop.size=pop.n,						#starting total population size
                                    start.full.pop.struct=c(),					#starting full pop structure - for if you want to switch off vacc
                                    vacc.t.multip=100*100,						#time to multiply out the vacc by chosen value - set to far beyond Tmax!!
                                    value.multip.vacc=1,						#value for 
                                    gamma=0.97)#,maxtest=75,step.size=ss)#, tol=0.1) ## add if make binary
  
  
  ## run matched on incidence 
  tmp1 <- iterateSIR(age.classes=age.classes, mort=mort,	
                     waifw=adj$waifw,
                     alpha=0,discrete.seas=c(),	
                     waning.maternal.par=0.5,waning.immunity=0,
                     vaccination=c(rep(0,11),vacc.cover,rep(0,116)),						
                     time.step=time.step,
                     fert= fert,alpha.births =0,
                     vacc.t.multip=vacc.t.multip,value.multip.vacc=value.multip.vacc,	#switch off vaccination from vacc.t.multip
                     burnin=(n.per.year*nyrs.burnin),Tmax=(n.per.year*nyrs), 
                     do.plot=F,start.pop=c(),start.pop.size=pop.n,start.full.pop.struct=c(),gamma=0.97)
  
  age.vals1 <- extractAge(tmp1$nt.store[tmp1$i.inds,],tmp1$age.classes)
  
  
  ## match age to previous
  adj2 <- tweakWAIFWtoMatchAge(target.age=mean(age.vals$avg.age[(nyrs.match*n.per.year):((nyrs.match+1)*n.per.year)]),                 
                                    test.multip=seq(0.1,5, length=10),
                                    age.classes=tmp$age.classes,  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                                    mort=mort,					#neglgible death till age 50 (could make more reaslitic)
                                    waifw=scaled.waifw.mat,					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                                    alpha=0,							#parameter governing sine wave seasonality in transmission	
                                    discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                                    waning.maternal.par=0.5,						#rate of waning of maternal immunity
                                    waning.immunity=0, 						#rate of waning of immunity in general
                                    vaccination=c(rep(0,11),vacc.cover,rep(0,116)),						
                                    time.step=time.step,							#time-step (1/4 of a month = 1 week)
                                    fert=fert,					#pattern of fertility over age
                                    alpha.births = 0,						#pattern of seasonality in births	
                                    tmatch=(nyrs.match*n.per.year):((nyrs.match+1)*n.per.year),burnin=500,Tmax=(nyrs.match+3)*n.per.year, 
                                    do.plot=FALSE,	
                                    start.pop=c(),							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                                    start.pop.size=pop.n,						#starting total population size
                                    start.full.pop.struct=c(),					#starting full pop structure - for if you want to switch off vacc
                                    vacc.t.multip=100*100,						#time to multiply out the vacc by chosen value - set to far beyond Tmax!!
                                    value.multip.vacc=1,						#value for 
                                    gamma=0.97)
  
  
  if(adj2$best==0.1) { new.range <- seq(0.01,0.1,length=20); ss <- 0.01}
  if(adj2$best==5) { new.range <- seq(5,20,length=20); ss <- 0.25}
  if(adj2$best!=5 & adj2$best!=0.1) {new.range <- seq(adj2$best-diff(adj2$test.multip)[1],adj2$best+diff(adj2$test.multip)[1],length=20); ss <- 0.01}
  first.best.age <- adj2$best
  
  print(c("ranges",first.best.incidence,first.best.age))
  
  ## match age to previous
  adj2 <- tweakWAIFWtoMatchAge(target.age=mean(age.vals$avg.age[(nyrs.match*n.per.year):((nyrs.match+1)*n.per.year)]),                 
#  adj2 <- tweakWAIFWtoMatchAgeBinary(target.age=mean(age.vals$avg.age[(nyrs.match*n.per.year):((nyrs.match+1)*n.per.year)]),                 
                               test.multip=new.range,#[c(1,20)],
                               age.classes=tmp$age.classes,  	#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                               mort=mort,					#neglgible death till age 50 (could make more reaslitic)
                               waifw=scaled.waifw.mat,					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                               alpha=0,							#parameter governing sine wave seasonality in transmission	
                               discrete.seas=c(),						#discrete seasonality from a TSIR fit - over-rides the previous	
                               waning.maternal.par=0.5,						#rate of waning of maternal immunity
                               waning.immunity=0, 						#rate of waning of immunity in general
                               vaccination=c(rep(0,11),vacc.cover,rep(0,116)),						
                               time.step=time.step,							#time-step (1/4 of a month = 1 week)
                               fert=fert,					#pattern of fertility over age
                               alpha.births = 0,						#pattern of seasonality in births	
                               tmatch=(nyrs.match*n.per.year):((nyrs.match+1)*n.per.year),burnin=500,Tmax=(nyrs.match+3)*n.per.year, 
                               do.plot=FALSE,	
                               start.pop=c(),							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start # 
                               start.pop.size=pop.n,						#starting total population size
                               start.full.pop.struct=c(),					#starting full pop structure - for if you want to switch off vacc
                               vacc.t.multip=100*100,						#time to multiply out the vacc by chosen value - set to far beyond Tmax!!
                               value.multip.vacc=1,						#value for 
                               gamma=0.97)#,maxtest=75,step.size=ss)

  
  ## run matched on incidence 
  tmp2 <- iterateSIR(age.classes=age.classes, mort=mort,	
                   waifw=adj2$waifw,
                   alpha=0,
                   waning.maternal.par=0.5,waning.immunity=0,
                   vaccination=c(rep(0,11),vacc.cover,rep(0,116)),						
                   time.step=time.step,
                   fert= fert,
                   vacc.t.multip=vacc.t.multip,value.multip.vacc=value.multip.vacc,	#switch off vaccination from vacc.t.multip
                   burnin=(n.per.year*nyrs.burnin),Tmax=(n.per.year*nyrs), 
                   do.plot=F,start.pop=c(),start.pop.size=pop.n,start.full.pop.struct=c(),gamma=0.97)

  age.vals2 <- extractAge(tmp2$nt.store[tmp2$i.inds,],tmp2$age.classes)


  
  ## get summaries and get ready
  
  summary1 <- extractStats(incidence=colSums(tmp$nt.store[tmp$i.inds,(vacc.t.multip):tmp$Tmax]),avg.age=age.vals$avg.age[(vacc.t.multip):tmp$Tmax],var.age=age.vals$var.age[(vacc.t.multip):tmp$Tmax]) 
  summary2 <- extractStats(incidence=colSums(tmp1$nt.store[tmp1$i.inds,(vacc.t.multip):tmp1$Tmax]),avg.age=age.vals1$avg.age[(vacc.t.multip):tmp$Tmax],var.age=age.vals1$var.age[(vacc.t.multip):tmp1$Tmax]) 
  summary3 <- extractStats(incidence=colSums(tmp2$nt.store[tmp2$i.inds,(vacc.t.multip):tmp2$Tmax]),avg.age=age.vals2$avg.age[(vacc.t.multip):tmp$Tmax],var.age=age.vals2$var.age[(vacc.t.multip):tmp2$Tmax]) 
  
  match.times <- (nyrs.match*n.per.year):((nyrs.match+1)*n.per.year)
  final.data <- data.frame(target.R0=rep(R0,3),type.waifw=c("flat","struct","struct"),type.match=c("none","incidence","age"),
                           tpeak=c(summary1$tpeak,summary2$tpeak,summary3$tpeak),
                           incidence.peak =c(max(colSums(tmp$nt.store[tmp$i.inds,(vacc.t.multip):tmp$Tmax])),max(colSums(tmp1$nt.store[tmp1$i.inds,(vacc.t.multip):tmp1$Tmax])),max(colSums(tmp2$nt.store[tmp1$i.inds,(vacc.t.multip):tmp1$Tmax]))),                           
                           incidence.pre.event=c(sum(tmp$nt.store[tmp$i.inds,match.times]),sum(tmp1$nt.store[tmp$i.inds,match.times]),sum(tmp2$nt.store[tmp$i.inds,match.times])),
                           incidence.at.end=c(sum(tmp$nt.store[tmp$i.inds,tmp$Tmax]),sum(tmp1$nt.store[tmp$i.inds,tmp1$Tmax]),sum(tmp2$nt.store[tmp$i.inds,tmp1$Tmax])),
			   cumulative.incidence = c(sum(tmp$nt.store[tmp$i.inds,summary1$tpeak:tmp$Tmax]),sum(tmp1$nt.store[tmp$i.inds,summary2$tpeak:tmp1$Tmax]),sum(tmp2$nt.store[tmp$i.inds,summary3$tpeak:tmp1$Tmax])),
                           age.pre.event=c(mean(age.vals$avg.age[match.times]),mean(age.vals1$avg.age[match.times]),mean(age.vals2$avg.age[match.times])),
                           age.at.end=c(age.vals$avg.age[tmp$Tmax],age.vals1$avg.age[tmp$Tmax],age.vals2$avg.age[tmp$Tmax]),
                           R0=c(getR0(scaled.flat.waifw,xx$stable.age*pop.n),getR0(adj$waifw,xx$stable.age*pop.n),getR0(adj2$waifw,xx$stable.age*pop.n)),
                           max.waifw=c(max(scaled.flat.waifw),max(adj$waifw),max(adj2$waifw)))
                      
  

   
  if (full.plot){
   par(mfrow=c(1,2), bty="l")
   image(age.classes/12,age.classes/12, adj$waifw, xlab="Age (years)",ylab="Age (years)")
   plot(tmp$age.classes/12,tmp$nt.store[tmp$i.inds,tmp$Tmax], type="l", xlim=c(0,10), xlab="Age (years)", ylab="Incidence",lwd=2)
   points(tmp$age.classes/12,tmp1$nt.store[tmp$i.inds,tmp$Tmax], type="l", col=4,lwd=2)
   points(tmp$age.classes/12,tmp2$nt.store[tmp$i.inds,tmp$Tmax], type="l", col="forestgreen",lwd=2)
  
  
  ## 
  par(mfrow=c(4,2), mar=c(5,4,3,2))
  plot((tmp$burnin:tmp$Tmax)/n.per.year,colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]), type="l", xlab="time", ylab="incidence", ylim=range(c(0,colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]))))
  abline(v=tmp$vacc.t.multip/n.per.year,lty=3,col="red")
  title("Flat waifw")
  
  plot(colSums(tmp$nt.store[tmp$s.inds,tmp$burnin:tmp$Tmax]), colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]),type="l",  ylab="incidence", xlab="susceptible")
  points(colSums(tmp$nt.store[tmp$s.inds,tmp$vacc.t.multip:tmp$Tmax]), colSums(tmp$nt.store[tmp$i.inds,tmp$vacc.t.multip:tmp$Tmax]),type="l", col="black")
  
  plot((tmp$burnin:tmp$Tmax)/n.per.year,colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]), type="l", xlab="time", ylab="incidence", ylim=range(c(0,colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]))), col="blue")
  abline(v=tmp$vacc.t.multip/n.per.year,lty=3,col="red")
  title("Diagonal waifw, incidence matched")

  plot(colSums(tmp1$nt.store[tmp$s.inds,tmp$burnin:tmp$Tmax]), colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]),type="l", ylab="incidence", xlab="susceptible", col="blue")
  points(colSums(tmp1$nt.store[tmp$s.inds,tmp$vacc.t.multip:tmp$Tmax]), colSums(tmp1$nt.store[tmp$i.inds,tmp$vacc.t.multip:tmp$Tmax]),type="l", col="blue")

  plot((tmp$burnin:tmp$Tmax)/n.per.year,colSums(tmp2$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]), type="l", xlab="time", ylab="incidence", ylim=range(c(0,colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]))), col="forestgreen")
  abline(v=tmp$vacc.t.multip/n.per.year,lty=3,col="red")
  title("Diagonal waifw, age matched")
  
  plot(colSums(tmp2$nt.store[tmp$s.inds,tmp$burnin:tmp$Tmax]), colSums(tmp2$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]),type="l", ylab="incidence", xlab="susceptible", col="forestgreen")
  points(colSums(tmp2$nt.store[tmp$s.inds,tmp$vacc.t.multip:tmp$Tmax]), colSums(tmp2$nt.store[tmp$i.inds,tmp$vacc.t.multip:tmp$Tmax]),type="l", col="forestgreen")
  
    
  plot((tmp$burnin:tmp$Tmax)/n.per.year,age.vals$avg.age[tmp$burnin:tmp$Tmax]/12, type="l", ylim=range(c(age.vals$avg.age[tmp$burnin:tmp$Tmax],age.vals1$avg.age[tmp$burnin:tmp$Tmax],age.vals2$avg.age[tmp$burnin:tmp$Tmax])/12), ylab="Avg age of infection", xlab="time")
  points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals1$avg.age[tmp$burnin:tmp$Tmax]/12, type="l", col="blue")
  points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals2$avg.age[tmp$burnin:tmp$Tmax]/12, type="l", col="forestgreen")
  abline(v=tmp$vacc.t.multip/n.per.year,lty=3,col="red")
  
  plot((tmp$burnin:tmp$Tmax)/n.per.year,age.vals$var.age[tmp$burnin:tmp$Tmax]/12, type="l", ylim=range(c(age.vals$var.age[tmp$burnin:tmp$Tmax],age.vals1$var.age[tmp$burnin:tmp$Tmax],age.vals2$var.age[tmp$burnin:tmp$Tmax])/12), ylab="Var age of infection", xlab="time")
  points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals1$var.age[tmp$burnin:tmp$Tmax]/12, type="l", col="blue")
  points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals2$var.age[tmp$burnin:tmp$Tmax]/12, type="l", col="forestgreen")
  abline(v=tmp$vacc.t.multip/n.per.year,lty=3,col="red")
  
  ### if you have the same R0 should you have the same # of susceptibles? definition is # new infections in a completely susceptible population for the avg infected individual
  ### could match infected instead? 

  
  par(mfrow=c(3,1), mar=c(5,4,3,2))
  plot((tmp$burnin:tmp$Tmax)/n.per.year,colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]), type="l", xlab="time", ylab="incidence", 
       ylim=range(c(0,colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]),
                      colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]),
                      colSums(tmp2$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]))))
  points((tmp$burnin:tmp$Tmax)/n.per.year,colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]), type="l",col="blue")
  points((tmp$burnin:tmp$Tmax)/n.per.year,colSums(tmp2$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]), type="l",col="forestgreen")
  abline(v=tmp$vacc.t.multip/n.per.year,lty=3,col="red")
  
  plot((tmp$burnin:tmp$Tmax)/n.per.year,age.vals$avg.age[tmp$burnin:tmp$Tmax]/12, type="l", ylim=range(c(age.vals$avg.age[tmp$burnin:tmp$Tmax],age.vals1$avg.age[tmp$burnin:tmp$Tmax],age.vals2$avg.age[tmp$burnin:tmp$Tmax])/12), ylab="Avg age of infection", xlab="time")
  points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals1$avg.age[tmp$burnin:tmp$Tmax]/12, type="l", col="blue")
  points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals2$avg.age[tmp$burnin:tmp$Tmax]/12, type="l", col="forestgreen")
  abline(v=tmp$vacc.t.multip/n.per.year,lty=3,col="red")
  
  plot((tmp$burnin:tmp$Tmax)/n.per.year,age.vals$var.age[tmp$burnin:tmp$Tmax]/12, type="l", ylim=range(c(age.vals$var.age[tmp$burnin:tmp$Tmax],age.vals1$var.age[tmp$burnin:tmp$Tmax],age.vals2$var.age[tmp$burnin:tmp$Tmax])/12), ylab="Var age of infection", xlab="time")
  points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals1$var.age[tmp$burnin:tmp$Tmax]/12, type="l", col="blue")
  points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals2$var.age[tmp$burnin:tmp$Tmax]/12, type="l", col="forestgreen")
  abline(v=tmp$vacc.t.multip/n.per.year,lty=3,col="red")

  
  ##match on AGE pushes WAIFW way down, but the age being higher places it in a higher part of transmission at the end. 
  par(mfrow=c(1,1))
  plot(tmp$age.classes/12,colSums(tmp1$Tmat[tmp$i.inds,tmp1$s.inds]),  xlab="Age", ylab="Transitions from S to I", col="blue", type="l", lwd=2)
  points(tmp$age.classes/12,colSums(tmp2$Tmat[tmp$i.inds,tmp1$s.inds]),type="l", col="forestgreen", lwd=2)
  points(tmp$age.classes/12,colSums(tmp$Tmat[tmp$i.inds,tmp1$s.inds]),type="l",  lwd=2)
  
  abline(v=c(age.vals1$avg.age[tmp$burnin]/12,age.vals2$avg.age[tmp$burnin]/12), col=c("blue","forestgreen"), lty=3)
  abline(v=c(age.vals1$avg.age[tmp$Tmax]/12,age.vals2$avg.age[tmp$Tmax]/12), col=c("blue","forestgreen"))
  
  }
  
return(final.data)
  
}





### wrapper model, also adjusted to be measles like parameters, aka, time-step =1/2, and 26 where there is 52
### This function allows you to pck an R0 and either add or remove vaccination for a particular waifw 
###    It will start by running dynamics for a flat waifw set to match the R0
###     Then will match the dynamics for the new waifw on either incidence or age to the time before the intervention (remove / introduce vaccination)

aFullComparisonSimple <- function(R0=30, remove.vacc= TRUE, 
                            waifw.test=get.smooth.WAIFW(age.class.boundries=c(1:120,seq(180,600,by=60))/12, mu=30,sig=1, gam=0.08, delta=0.01), #blobby
			    index=1,		## index i-v captures the five different waifws being deployed, as an index to track down susceptibles and incidence
                            full.plot=FALSE, do.unity=FALSE){
  
  n.per.year <- 26
  time.step <- 1/2
  
  #starting parameters
  nyrs <- 100
  nyrs.burnin <- 70
  nyrs.match <- 70            #take this +1 to match 
  vacc.t.multip <- 80*n.per.year 
  pop.n <- 500000
  ## remove vacc
  if (remove.vacc) { 
    vacc.cover <- 0.7
    value.multip.vacc <- 0 
    value.add.vacc <- 0
    } else {
      ## introduce vacc
      vacc.cover <-0   ## since defined as a multiplier, original must be >0
      value.multip.vacc <- 1
      value.add.vacc <- 0.7
  }
  age.classes <- c(1:120,seq(180,600,by=60))
  mort <- c(rep(0,127),1)
  fert <-  c(rep(0,127),1) ##note that if lambda is not 1; vast slippage in average age
  
  #waning.maternal <- 0.99 #very rapid TODO CHECK
  waning.maternal <- 0.5 
 
  ## get a year of times to do the matching on age / incidence
  match.times <- (nyrs.match*n.per.year):((nyrs.match+1)*n.per.year)
  
  ## tweak this outside of this function to adjust the rate of increase to get something close to equilbirium (corresponds to lambda=1)
  ## here, using this to get reasonable starting values
  xx <- findStableStruct(age.classes=age.classes, 
                         mort=mort,fert =fert,time.step = time.step)
  #xx$lambda

  ## tune the flat waifw to get the target R0 
  scaled.flat.waifw <- scale.Waifw(R0=R0, DFE=xx$stable.age*pop.n, waifw=matrix(0.0002,length(age.classes),length(age.classes)))
  


  ## 1. Flat Waifw
  tmp <- iterateSIR(age.classes=age.classes, mort=mort,	
                    waifw=scaled.flat.waifw,
                    alpha=0,
                    waning.maternal.par=waning.maternal,waning.immunity=0,
                    vaccination=c(rep(0,11),vacc.cover,rep(0,116)),						
                    time.step=time.step,
                    fert= fert,
                    vacc.t.multip=100*1000, 			#CHANGE should be far in future, no? 
		    value.multip.vacc=value.multip.vacc,	#switch off (or on) vaccination from vacc.t.multip
		    value.add.vacc=value.add.vacc, 			
                    burnin=(n.per.year*nyrs.burnin),
		    Tmax=(n.per.year*nyrs), 
                    do.plot=F,start.pop=c(),start.pop.size=pop.n,
		    start.full.pop.struct=c(),gamma=0.97, do.unity=do.unity)
  
  age.vals <- extractAge(tmp$nt.store[tmp$i.inds,],tmp$age.classes)
  


  
  ## 2. Target waifw - scaled initially using the R0 trick 
  scaled.waifw.mat <- scale.Waifw(R0=R0, DFE=xx$stable.age*pop.n, waifw=waifw.test)
  
  ##
  #print("TARGET")
  #print(sum(tmp$nt.store[tmp$i.inds,match.times]))
  #print(mean(age.vals$avg.age[match.times]))
  ##  

  adj <- tweakWAIFWtoMatchBoth(target.incidence=sum(tmp$nt.store[tmp$i.inds,match.times]),
                               target.age=mean(age.vals$avg.age[match.times]),   
                               test.multip=c(seq(0.1,3,by=0.05),seq(3.5,31,length=20),39), 			#the set of amplitudes to test 
                               age.classes=tmp$age.classes,  				#up to 5 years in months, up to 10 in years, and up to 50 in 5 years 
                               mort=mort,						#neglgible death till age 50 (could make more reaslitic)
                               waifw=scaled.waifw.mat,					#transmission over age - set to the same value everywhere equal to transmission rate if want no structure!
                               alpha=0,							#parameter governing sine wave seasonality in transmission	
                               discrete.seas=c(),					#discrete seasonality from a TSIR fit - over-rides the previous	
                               waning.maternal.par=waning.maternal,			#rate of waning of maternal immunity
                               waning.immunity=0, 					#rate of waning of immunity in general
                               vaccination=c(rep(0,11),vacc.cover,rep(0,116)),						
                               time.step=time.step,					#time-step 
                               fert=fert,						#pattern of fertility over age
                               alpha.births = 0,						#pattern of seasonality in births	
                               tmatch=(nyrs.match*n.per.year):((nyrs.match+1)*n.per.year),
			       burnin=500,Tmax=(nyrs.match+3)*n.per.year, 
                               do.plot=FALSE,	
                               start.pop=c(),							#starting population structure and # - could impose, otherwise takes the equilinrium structure with imposed start 
                               start.pop.size=pop.n,						#starting total population size
                               start.full.pop.struct=c(),					#starting full pop structure 
                               vacc.t.multip=100*100,						#time to multiply out the vacc by chosen value - set to far beyond Tmax!! since matching on initial hear
                               value.multip.vacc=1,
			       value.add.vacc=0, 						
                               gamma=0.97)
  
  
  ## run matched on incidence 
  tmp1 <- iterateSIR(age.classes=age.classes, mort=mort,	
                     waifw=adj$waifw,
                     alpha=0,discrete.seas=c(),	
                     waning.maternal.par=waning.maternal,waning.immunity=0,
                     vaccination=c(rep(0,11),vacc.cover,rep(0,116)),						
                     time.step=time.step,
                     fert= fert,alpha.births =0,
                     vacc.t.multip=vacc.t.multip,value.multip.vacc=value.multip.vacc,	#switch off vaccination from vacc.t.multip
		     value.add.vacc=value.add.vacc, 
                     burnin=(n.per.year*nyrs.burnin),Tmax=(n.per.year*nyrs), 
                     do.plot=F,start.pop=c(),start.pop.size=pop.n,start.full.pop.struct=c(),gamma=0.97, do.unity=do.unity)
  
  age.vals1 <- extractAge(tmp1$nt.store[tmp1$i.inds,],tmp1$age.classes)
  #print(sum(tmp1$nt.store[tmp$i.inds,match.times]))
  
  
  ## run matched on age
  tmp2 <- iterateSIR(age.classes=age.classes, mort=mort,	
                     waifw=adj$waifw.age,
                     alpha=0,
                     waning.maternal.par=waning.maternal,waning.immunity=0,
                     vaccination=c(rep(0,11),vacc.cover,rep(0,116)),						
                     time.step=time.step,
                     fert= fert,
                     vacc.t.multip=vacc.t.multip,value.multip.vacc=value.multip.vacc,	#switch off vaccination from vacc.t.multip
                      value.add.vacc=value.add.vacc, 
                     burnin=(n.per.year*nyrs.burnin),Tmax=(n.per.year*nyrs), 
                     do.plot=F,start.pop=c(),start.pop.size=pop.n,start.full.pop.struct=c(),gamma=0.97, do.unity=do.unity)
  
  age.vals2 <- extractAge(tmp2$nt.store[tmp2$i.inds,],tmp2$age.classes)
  
  
  
  ## get summaries and get ready
  
  summary1 <- extractStats(incidence=colSums(tmp$nt.store[tmp$i.inds,(vacc.t.multip):tmp$Tmax]),avg.age=age.vals$avg.age[(vacc.t.multip):tmp$Tmax],var.age=age.vals$var.age[(vacc.t.multip):tmp$Tmax]) 
  summary2 <- extractStats(incidence=colSums(tmp1$nt.store[tmp1$i.inds,(vacc.t.multip):tmp1$Tmax]),avg.age=age.vals1$avg.age[(vacc.t.multip):tmp$Tmax],var.age=age.vals1$var.age[(vacc.t.multip):tmp1$Tmax]) 
  summary3 <- extractStats(incidence=colSums(tmp2$nt.store[tmp2$i.inds,(vacc.t.multip):tmp2$Tmax]),avg.age=age.vals2$avg.age[(vacc.t.multip):tmp$Tmax],var.age=age.vals2$var.age[(vacc.t.multip):tmp2$Tmax]) 
  
  print(summary1$tpeak)
  
  final.data <- data.frame(target.R0=rep(R0,3),type.waifw=c("flat","struct","struct"),type.match=c("none","incidence","age"),
                           tpeak=c(summary1$tpeak,summary2$tpeak,summary3$tpeak),
                           incidence.peak =c(max(colSums(tmp$nt.store[tmp$i.inds,(vacc.t.multip):tmp$Tmax])),max(colSums(tmp1$nt.store[tmp1$i.inds,(vacc.t.multip):tmp1$Tmax])),max(colSums(tmp2$nt.store[tmp1$i.inds,(vacc.t.multip):tmp1$Tmax]))),                           
                           incidence.pre.event=c(sum(tmp$nt.store[tmp$i.inds,match.times]),sum(tmp1$nt.store[tmp1$i.inds,match.times]),sum(tmp2$nt.store[tmp$i.inds,match.times])),
                           incidence.at.end=c(sum(tmp$nt.store[tmp$i.inds,tmp$Tmax]),sum(tmp1$nt.store[tmp$i.inds,tmp1$Tmax]),sum(tmp2$nt.store[tmp$i.inds,tmp1$Tmax])),
  			   cumulative.incidence = c(sum(tmp$nt.store[tmp$i.inds,summary1$tpeak:tmp$Tmax]),sum(tmp1$nt.store[tmp$i.inds,summary2$tpeak:tmp1$Tmax]),sum(tmp2$nt.store[tmp$i.inds,summary3$tpeak:tmp1$Tmax])),
                           age.pre.event=c(mean(age.vals$avg.age[match.times]),mean(age.vals1$avg.age[match.times]),mean(age.vals2$avg.age[match.times])),
                           age.after.five.years=c(age.vals$avg.age[vacc.t.multip+5*n.per.year],age.vals1$avg.age[vacc.t.multip+5*n.per.year],age.vals2$avg.age[vacc.t.multip+5*n.per.year]),
                           age.at.end=c(age.vals$avg.age[tmp$Tmax],age.vals1$avg.age[tmp$Tmax],age.vals2$avg.age[tmp$Tmax]),
                           susceptible.pre.event=c(sum(tmp$nt.store[tmp$s.inds,match.times]),sum(tmp1$nt.store[tmp1$s.inds,match.times]),sum(tmp2$nt.store[tmp$s.inds,match.times])),
                           susceptible.at.end=c(sum(tmp$nt.store[tmp$s.inds,tmp$Tmax]),sum(tmp1$nt.store[tmp$s.inds,tmp1$Tmax]),sum(tmp2$nt.store[tmp$s.inds,tmp1$Tmax])),
                           R0=c(getR0(scaled.flat.waifw,xx$stable.age*pop.n),getR0(adj$waifw,xx$stable.age*pop.n),getR0(adj$waifw.age,xx$stable.age*pop.n)),
                           max.waifw=c(max(scaled.flat.waifw),max(adj$waifw),max(adj$waifw.age)))
  
  ## note that age.pre.event is a average while age.age.end and age.after.five.years are both at a single time point
  
  
  if (full.plot){
    par(mfrow=c(1,2), bty="l")
    image(age.classes/12,age.classes/12, adj$waifw, xlab="Age (years)",ylab="Age (years)")
    plot(tmp$age.classes/12,tmp$nt.store[tmp$i.inds,tmp$Tmax], type="l", xlim=c(0,10), xlab="Age (years)", ylab="Incidence",lwd=2)
    points(tmp$age.classes/12,tmp1$nt.store[tmp$i.inds,tmp$Tmax], type="l", col=4,lwd=2)
    points(tmp$age.classes/12,tmp2$nt.store[tmp$i.inds,tmp$Tmax], type="l", col="forestgreen",lwd=2)
    
    
    ## 
    par(mfrow=c(4,2), mar=c(5,4,3,2))
    plot((tmp$burnin:tmp$Tmax)/n.per.year,colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]), type="l", xlab="time", ylab="incidence", ylim=range(c(0,colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]))))
    abline(v=tmp1$vacc.t.multip/n.per.year,lty=3,col="red")
    title("Flat waifw")
    
    #phase plane uniform
    plot(colSums(tmp$nt.store[tmp$s.inds,tmp$burnin:tmp$Tmax]), colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]),type="l",  ylab="incidence", xlab="susceptible")
    points(colSums(tmp$nt.store[tmp$s.inds,tmp1$vacc.t.multip:tmp$Tmax]), colSums(tmp$nt.store[tmp$i.inds,tmp1$vacc.t.multip:tmp$Tmax]),type="l", col="black")
    
    plot((tmp$burnin:tmp$Tmax)/n.per.year,colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]), type="l", xlab="time", ylab="incidence", 
         ylim=range(c(0,colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]))), col="blue")
    abline(v=tmp1$vacc.t.multip/n.per.year,lty=3,col="red")
    title("Diagonal waifw, incidence matched")
    
    #phase plane uniform
    plot(colSums(tmp1$nt.store[tmp$s.inds,tmp$burnin:tmp$Tmax]), colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]),type="l", ylab="incidence", xlab="susceptible", col="blue")
    points(colSums(tmp1$nt.store[tmp$s.inds,tmp1$vacc.t.multip:tmp$Tmax]), colSums(tmp1$nt.store[tmp$i.inds,tmp1$vacc.t.multip:tmp$Tmax]),type="l", col="blue")
    
    plot((tmp$burnin:tmp$Tmax)/n.per.year,colSums(tmp2$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]), type="l", xlab="time", ylab="incidence", 
		ylim=range(c(0,colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]))), col="forestgreen")
    abline(v=tmp1$vacc.t.multip/n.per.year,lty=3,col="red")
    title("Diagonal waifw, age matched")
    
    #phase plane match on cidence
    plot(colSums(tmp2$nt.store[tmp$s.inds,tmp1$burnin:tmp$Tmax]), colSums(tmp2$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]),type="l", ylab="incidence", xlab="susceptible", col="forestgreen")
    points(colSums(tmp2$nt.store[tmp$s.inds,tmp1$vacc.t.multip:tmp$Tmax]), colSums(tmp2$nt.store[tmp$i.inds,tmp1$vacc.t.multip:tmp$Tmax]),type="l", col="forestgreen")
    
    
    plot((tmp$burnin:tmp$Tmax)/n.per.year,age.vals$avg.age[tmp$burnin:tmp$Tmax]/12, type="l", 
         ylim=range(c(age.vals$avg.age[tmp$burnin:tmp$Tmax],age.vals1$avg.age[tmp$burnin:tmp$Tmax],age.vals2$avg.age[tmp$burnin:tmp$Tmax])/12,na.rm=TRUE), ylab="Avg age of infection", xlab="time")
    points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals1$avg.age[tmp$burnin:tmp$Tmax]/12, type="l", col="blue")
    points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals2$avg.age[tmp$burnin:tmp$Tmax]/12, type="l", col="forestgreen")
    abline(v=tmp1$vacc.t.multip/n.per.year,lty=3,col="red")
    
    #phase plane match on age
    plot((tmp$burnin:tmp$Tmax)/n.per.year,age.vals$var.age[tmp$burnin:tmp$Tmax]/12, type="l", 
         ylim=range(c(age.vals$var.age[tmp$burnin:tmp$Tmax],age.vals1$var.age[tmp$burnin:tmp$Tmax],age.vals2$var.age[tmp$burnin:tmp$Tmax])/12,na.rm=TRUE), ylab="Var age of infection", xlab="time")
    points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals1$var.age[tmp$burnin:tmp$Tmax]/12, type="l", col="blue")
    points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals2$var.age[tmp$burnin:tmp$Tmax]/12, type="l", col="forestgreen")
    abline(v=tmp1$vacc.t.multip/n.per.year,lty=3,col="red")
    
    ### if you have the same R0 should you have the same # of susceptibles? definition is # new infections in a completely susceptible population for the avg infected individual
    ### could match infected instead? 
    
    
    par(mfrow=c(3,1), mar=c(5,4,3,2))
    plot((tmp$burnin:tmp$Tmax)/n.per.year,colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]), type="l", xlab="time", ylab="incidence", 
         ylim=range(c(0,colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]),
                      colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]),
                      colSums(tmp2$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax])),na.rm=TRUE))
    points((tmp$burnin:tmp$Tmax)/n.per.year,colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]), type="l",col="blue")
    points((tmp$burnin:tmp$Tmax)/n.per.year,colSums(tmp2$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]), type="l",col="forestgreen")
    abline(v=tmp1$vacc.t.multip/n.per.year,lty=3,col="red")
    
    plot((tmp$burnin:tmp$Tmax)/n.per.year,age.vals$avg.age[tmp$burnin:tmp$Tmax]/12, type="l", 
         ylim=range(c(age.vals$avg.age[tmp$burnin:tmp$Tmax],age.vals1$avg.age[tmp$burnin:tmp$Tmax],age.vals2$avg.age[tmp$burnin:tmp$Tmax])/12,na.rm=TRUE),
         ylab="Avg age of infection", xlab="time")
    points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals1$avg.age[tmp$burnin:tmp$Tmax]/12, type="l", col="blue")
    points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals2$avg.age[tmp$burnin:tmp$Tmax]/12, type="l", col="forestgreen")
    abline(v=tmp1$vacc.t.multip/n.per.year,lty=3,col="red")
    
    plot((tmp$burnin:tmp$Tmax)/n.per.year,age.vals$var.age[tmp$burnin:tmp$Tmax]/12, type="l", 
         ylim=range(c(age.vals$var.age[tmp$burnin:tmp$Tmax],age.vals1$var.age[tmp$burnin:tmp$Tmax],age.vals2$var.age[tmp$burnin:tmp$Tmax])/12,na.rm=TRUE), ylab="Var age of infection", xlab="time")
    points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals1$var.age[tmp$burnin:tmp$Tmax]/12, type="l", col="blue")
    points((tmp$burnin:tmp$Tmax)/n.per.year,age.vals2$var.age[tmp$burnin:tmp$Tmax]/12, type="l", col="forestgreen")
    abline(v=tmp1$vacc.t.multip/n.per.year,lty=3,col="red")
    
    
    ##match on AGE pushes WAIFW way down, but the age being higher places it in a higher part of transmission at the end. 
    par(mfrow=c(1,1))
    plot(tmp$age.classes/12,colSums(tmp1$Tmat[tmp$i.inds,tmp1$s.inds]),  xlab="Age", ylab="Transitions from S to I", col="blue", type="l", lwd=2)
    points(tmp$age.classes/12,colSums(tmp2$Tmat[tmp$i.inds,tmp1$s.inds]),type="l", col="forestgreen", lwd=2)
    points(tmp$age.classes/12,colSums(tmp$Tmat[tmp$i.inds,tmp1$s.inds]),type="l",  lwd=2)
    
    abline(v=c(age.vals1$avg.age[tmp$burnin]/12,age.vals2$avg.age[tmp$burnin]/12), col=c("blue","forestgreen"), lty=3)
    abline(v=c(age.vals1$avg.age[tmp$Tmax]/12,age.vals2$avg.age[tmp$Tmax]/12), col=c("blue","forestgreen"))
    
  }
  
   ## this seems to be the only place that burnin crops up
   if (index==1 & R0==10) append <- FALSE else append <- TRUE
   incidences <- cbind(rep(index,3),rep(R0,3),
		rbind(colSums(tmp$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]),colSums(tmp1$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax]),colSums(tmp2$nt.store[tmp$i.inds,tmp$burnin:tmp$Tmax])))
   susceptibles <- cbind(rep(index,3),rep(R0,3),
		rbind(colSums(tmp$nt.store[tmp$s.inds,tmp$burnin:tmp$Tmax]),colSums(tmp1$nt.store[tmp$s.inds,tmp$burnin:tmp$Tmax]),colSums(tmp2$nt.store[tmp$s.inds,tmp$burnin:tmp$Tmax])))
   unity.beta <- cbind(rep(index,3),rep(R0,3),
		rbind(tmp$store.beta.unity[tmp$burnin:tmp$Tmax],tmp1$store.beta.unity[tmp$burnin:tmp$Tmax],tmp2$store.beta.unity[tmp$burnin:tmp$Tmax]))
   all.age.vals <- cbind(rep(index,3),rep(R0,3),
			rbind(age.vals$avg.age[tmp$burnin:tmp$Tmax],age.vals1$avg.age[tmp$burnin:tmp$Tmax],age.vals2$avg.age[tmp$burnin:tmp$Tmax]))
   all.var.age.vals <- cbind(rep(index,3),rep(R0,3),
			rbind(age.vals$var.age[tmp$burnin:tmp$Tmax],age.vals1$var.age[tmp$burnin:tmp$Tmax],age.vals2$var.age[tmp$burnin:tmp$Tmax]))

   index.multip <- cbind(rep(index,2),c(R0,R0),c(adj$best,adj$best.age))

	#print(dim(unity.beta))
	#print(unity.beta[,1:10])
	
   write.table(incidences, file=paste("/Users/cmetcalf/Documents/temp/incidence",remove.vacc,".csv", sep=""),append=append,  col.names = !append,sep=",") 
   write.table(susceptibles, file=paste("/Users/cmetcalf/Documents/temp/susceptibles",remove.vacc,".csv", sep=""),append=append,  col.names = !append, sep=",")
   if (do.unity) write.table(unity.beta,file=paste("/Users/cmetcalf/Documents/temp/unity.beta",remove.vacc,".csv", sep=""),  col.names = !append,append=append, sep=",")
   write.table(all.age.vals, file=paste("/Users/cmetcalf/Documents/temp/age.vals",remove.vacc,".csv", sep=""),append=append,  col.names = !append, sep=",")
   write.table(all.var.age.vals, file=paste("/Users/cmetcalf/Documents/temp/var.age.vals",remove.vacc,".csv", sep=""),append=append,  col.names = !append, sep=",")

   write.table(index.multip, file=paste("/Users/cmetcalf/Documents/temp/adj.vals",remove.vacc,".csv", sep=""),append=append,  col.names = !append, sep=",")


  return(final.data)
  
}


### CURRENT RESULTS
checkOutResults <- function(){

   remove.vacc <- FALSE # TRUE
   R0test <- seq(10,24,by=2)
   chs <- 1

   incidences <- read.csv(file=paste("/Users/cmetcalf/Documents/temp/incidence",remove.vacc,".csv", sep=""), header=T, row.names=NULL) 
   susceptibles <- read.csv(file=paste("/Users/cmetcalf/Documents/temp/susceptibles",remove.vacc,".csv", sep=""), header=T, row.names=NULL)
   store.beta.unity <- read.csv(file=paste("/Users/cmetcalf/Documents/temp/unity.beta",remove.vacc,".csv", sep=""), header=T, row.names=NULL)
   avg.age <- read.csv(file=paste("/Users/cmetcalf/Documents/temp/age.vals",remove.vacc,".csv", sep=""), header=T, row.names=NULL)
   var.age <- read.csv(file=paste("/Users/cmetcalf/Documents/temp/var.age.vals",remove.vacc,".csv", sep=""), header=T, row.names=NULL)

   index.waifw <- incidences[,2]; name.waifw <- c("Flat", "Peak Age 5", "Peak Age 10", "Diagonal-ish","POLYMOD")
   index.R0 <- incidences[,3]; index.match <- rep(1:3,length(index.waifw)/3) ## first is just running out time-series. 

   chs.times <- 4:length(susceptibles[1,]) #burnin has already been removed; just remove first columns containing R0, etc
   vacc.time <- 10*26 #starts 10 years after start 
   chs.times.after.vacc <- chs.times[chs.times>vacc.time]

  cols <- c(1,4,2) # blue is matched on incidence, red is matched on age; the first should be ignored - is just the uniform run for matching. 
  lwds <- c(2,1,1)	

   peak.incidence <- incidence.at.end <- incidence.at.start<-  array(dim=c(3,5,8))
   if (remove.vacc) g <- function(x) max(x) else g <- function(x) min(x)
      for (chs.R0 in 1:8) { 
       chs.row.baseline <- c(1:length(index.waifw))[index.waifw == 1 & index.R0==R0test[chs.R0]]
	peak.incidence[1,1,chs.R0] <- which(as.numeric(incidences[chs.row.baseline[2], chs.times.after.vacc])==g(as.numeric(incidences[chs.row.baseline[2], chs.times.after.vacc])))[1]
        incidence.at.end[1,1,chs.R0] <- incidences[chs.row.baseline[2], ncol(incidences)]
        incidence.at.start[1,1,chs.R0] <- incidences[chs.row.baseline[2], 4]
       for (k in 1:5) { 
       chs.row <- c(1:length(index.waifw))[index.waifw == k & index.R0==R0test[chs.R0]]	
	peak.incidence[2,k,chs.R0] <- which(as.numeric(incidences[chs.row[2], chs.times.after.vacc])==g(as.numeric(incidences[chs.row[2], chs.times.after.vacc])))[1]
 	peak.incidence[3,k,chs.R0] <- which(as.numeric(incidences[chs.row[3], chs.times.after.vacc])==g(as.numeric(incidences[chs.row[3], chs.times.after.vacc])))[1]
        incidence.at.end[2,k,chs.R0] <- incidences[chs.row[2], ncol(incidences)]
        incidence.at.end[3,k,chs.R0] <- incidences[chs.row[3], ncol(incidences)]
        incidence.at.start[2,k,chs.R0] <- incidences[chs.row[2], 4]
        incidence.at.start[3,k,chs.R0] <- incidences[chs.row[3], 4]
	}}



     par(mfrow=c(3,5))

      for (chs.R0 in 2:4) { 

       chs.row.baseline <- c(1:length(index.waifw))[index.waifw == 1 & index.R0==R0test[chs.R0]]	
       
     for (k in 1:5) { 
       chs.row <- c(1:length(index.waifw))[index.waifw == k & index.R0==R0test[chs.R0]]	
       plot(as.numeric(susceptibles[chs.row.baseline[2], chs.times]), as.numeric(incidences[chs.row.baseline[2], chs.times]), 
		type="l", xlab="susceptible", ylab="infected", col=1,lwd=2,
		ylim=range(incidences[chs.row,chs.times]), xlim=range(susceptibles[chs.row,chs.times]),main=name.waifw[index.waifw[chs.row[1]]])

     for (j in 2:length(chs.row)) points(as.numeric(susceptibles[chs.row[j], chs.times]), as.numeric(incidences[chs.row[j], chs.times]), 
		type="l",lwd=lwds[index.match[chs.row[j]]],col=cols[index.match[chs.row[j]]])
      }
  
      }


  par(mfrow=c(3,5))


      for (chs.R0 in 2:4) { 

       chs.row.baseline <- c(1:length(index.waifw))[index.waifw == 1 & index.R0==R0test[chs.R0]]	

     for (k in 1:5) { 
       chs.row <- c(1:length(index.waifw))[index.waifw == k & index.R0==R0test[chs.R0]]	
       plot(as.numeric(store.beta.unity[chs.row.baseline[2], chs.times]), type="l", xlab="", ylab=expression(hat(beta)), col=1,lwd=2,
		#ylim=range(store.beta.unity[chs.row,chs.times]),main=name.waifw[index.waifw[chs.row[j]]])
		ylim=range(store.beta.unity[index.R0==R0test[chs.R0],chs.times],na.rm=TRUE),main=name.waifw[index.waifw[chs.row[1]]])


     for (j in 2:length(chs.row)) points(as.numeric(store.beta.unity[chs.row[j], chs.times]),type="l",lwd=lwds[index.match[chs.row[j]]],col=cols[index.match[chs.row[j]]])
	abline(v=vacc.time,col="grey")

      }
  
      }



  par(mfrow=c(3,5))


      for (chs.R0 in 2:4) { 
     # for (chs.R0 in c(2,5,8)) { 

       chs.row.baseline <- c(1:length(index.waifw))[index.waifw == 1 & index.R0==R0test[chs.R0]]	

     for (k in 1:5) { 
       chs.row <- c(1:length(index.waifw))[index.waifw == k & index.R0==R0test[chs.R0]]	
       plot(as.numeric(avg.age[chs.row.baseline[2], chs.times])/12, type="l", xlab="", ylab="Average Age", col=1,lwd=2,
		ylim=range(c(avg.age[chs.row,chs.times],as.numeric(avg.age[chs.row.baseline[2], chs.times])),na.rm=TRUE)/12,main=name.waifw[index.waifw[chs.row[j]]])


     for (j in 2:length(chs.row)) points(as.numeric(avg.age[chs.row[j], chs.times])/12,type="l",lwd=lwds[index.match[chs.row[j]]],col=cols[index.match[chs.row[j]]])
	abline(v= vacc.time, col="grey")

      }
  
      }




  par(mfrow=c(3,5))


      for (chs.R0 in 2:4) { 
     # for (chs.R0 in c(2,5,8)) { 

       chs.row.baseline <- c(1:length(index.waifw))[index.waifw == 1 & index.R0==R0test[chs.R0]]	

     for (k in 1:5) { 
       chs.row <- c(1:length(index.waifw))[index.waifw == k & index.R0==R0test[chs.R0]]	
       plot(as.numeric(var.age[chs.row.baseline[2], chs.times])/12, type="l", xlab="", ylab="Variance Age", col=1,lwd=2,log="y",
		ylim=range(c(var.age[chs.row,chs.times],as.numeric(var.age[chs.row.baseline[2], chs.times])),na.rm=TRUE)/12,main=name.waifw[index.waifw[chs.row[j]]])


     for (j in 2:length(chs.row)) points(as.numeric(var.age[chs.row[j], chs.times])/12,type="l",lwd=lwds[index.match[chs.row[j]]],col=cols[index.match[chs.row[j]]])
	abline(v= vacc.time, col="grey")

      }
  
      }



matplot(t(incidences[incidences[,3]==14 & incidences[,1]==2,4:ncol(incidences)]), type="l", col=cols, lty=1)
matplot(t(susceptibles[incidences[,3]==14 & incidences[,1]==2,4:ncol(incidences)]),t(incidences[incidences[,3]==14 & incidences[,1]==2,4:ncol(incidences)]), type="l", col=cols, lty=1)


      ## get some benchmarks of differences to probe a bit 
      max.difference.self.incidence  <- max.difference.incidence <- start.difference.incidence <- end.difference.incidence <- matrix(NA,5,length(R0test))
      max.difference.self.age <- max.difference.age <- start.difference.age <- end.difference.age <- matrix(NA,5,length(R0test))
      organize.age.inc.end  <- organize.age.age.end <- organize.age.inc <- organize.age.age <- matrix(NA,5,length(R0test))
      organize.var.age.inc <- organize.var.age.age <- organize.var.age.inc.end <- organize.var.age.age.end <- matrix(NA,5,length(R0test))
      start.beta <- cumsum.beta.1yr  <- difference.self.incidence.1yr <- difference.self.age.1yr <- matrix(NA,5,length(R0test))

      for (chs.R0 in 1:length(R0test)) { 
	     chs.row.baseline <- c(1:length(index.waifw))[index.waifw == 1 & index.R0==R0test[chs.R0]]	

  	   for (k in 1:5) { 
   		 chs.row <- c(1:length(index.waifw))[index.waifw == k & index.R0==R0test[chs.R0]]	
   
		incidence.matched.differences <- as.numeric(store.beta.unity[chs.row.baseline[2], chs.times])- as.numeric(store.beta.unity[chs.row[2], chs.times])
		age.matched.differences <- as.numeric(store.beta.unity[chs.row.baseline[2], chs.times])- as.numeric(store.beta.unity[chs.row[3], chs.times])
		incidence.matched.differences.self <-  as.numeric(store.beta.unity[chs.row[2], chs.times[1]])- as.numeric(store.beta.unity[chs.row[2], chs.times])
		age.matched.differences.self <- as.numeric(store.beta.unity[chs.row[2], chs.times[1]])- as.numeric(store.beta.unity[chs.row[3], chs.times])

		max.difference.incidence[k,chs.R0] <- median(incidence.matched.differences[which(abs(incidence.matched.differences)==max(abs(incidence.matched.differences)))],na.rm=TRUE)
		max.difference.age[k,chs.R0] <- median(age.matched.differences[which(abs(age.matched.differences)==max(abs(age.matched.differences)))],na.rm=TRUE)
		end.difference.incidence[k,chs.R0] <- incidence.matched.differences[length(incidence.matched.differences)]
		end.difference.age[k,chs.R0] <- age.matched.differences[length(age.matched.differences)]
		start.difference.incidence[k,chs.R0] <- incidence.matched.differences[1]
		start.difference.age[k,chs.R0] <- age.matched.differences[1]

		max.difference.self.incidence[k,chs.R0] <- -median(incidence.matched.differences.self[which(abs(incidence.matched.differences.self)==max(abs(incidence.matched.differences.self)))],na.rm=TRUE)
		max.difference.self.age[k,chs.R0] <- -median(age.matched.differences.self[which(abs(age.matched.differences.self)==max(abs(age.matched.differences.self)))],na.rm=TRUE)

		organize.age.inc[k,chs.R0] <- as.numeric(avg.age[chs.row[2], chs.times[1]])
		organize.age.age[k,chs.R0] <- as.numeric(avg.age[chs.row[3], chs.times[1]])

		organize.age.inc.end[k,chs.R0] <- as.numeric(avg.age[chs.row[2], chs.times[length(chs.times)]])
		organize.age.age.end[k,chs.R0] <- as.numeric(avg.age[chs.row[3], chs.times[length(chs.times)]])

		organize.var.age.inc[k,chs.R0] <- as.numeric(var.age[chs.row[2], chs.times[1]])
		organize.var.age.age[k,chs.R0] <- as.numeric(var.age[chs.row[3], chs.times[1]])

		organize.var.age.inc.end[k,chs.R0] <- as.numeric(var.age[chs.row[2], chs.times[length(chs.times)]])
		organize.var.age.age.end[k,chs.R0] <- as.numeric(var.age[chs.row[3], chs.times[length(chs.times)]])


		difference.self.incidence.1yr[k,chs.R0] <- -incidence.matched.differences.self[260+26] #10*26 years is when vaccination introduced, 26 is 1 year later
		difference.self.age.1yr[k,chs.R0] <- -age.matched.differences.self[260+26]

		cumsum.beta.1yr[k,chs.R0] <- cumsum(as.numeric(store.beta.unity[chs.row[2], chs.times]))[260+26]
		start.beta[k,chs.R0] <- cumsum(as.numeric(store.beta.unity[chs.row[2], chs.times]))[1]


	}}


	cols <- c(1,4,3,5,2)

	par(mfrow=c(1,2), mar=c(4,5,2,1), bty="l")
	#matplot(R0test, t(start.difference.incidence), ylab=expression(hat(beta[f])*" - "*hat(beta[a])), type="l", col=cols,lty=1, xlab=expression(R[0])); 
	#title("At start"); abline(h=0)
	legend("bottomright",legend= name.waifw, col=cols,lty=1,bty="n")
	#matplot(R0test, t(end.difference.incidence), ylab=expression(hat(beta)[f]*" - "*hat(beta[a])), type="l", col=cols,lty=1, xlab=expression(R[0])); title("At end"); abline(h=0)
	#matplot(R0test, t(max.difference.incidence), ylab=expression(hat(beta)[f]*" - "*hat(beta[a])), type="l", col=cols,lty=1, xlab=expression(R[0])); title("Largest (absolute)"); abline(h=0)
matplot(R0test, t(max.difference.self.incidence), ylab=expression(hat(beta)[a]*" - "*hat(beta[a1])), type="l", col=cols,lty=1, xlab=expression(R[0])); 
legend("topleft",legend= name.waifw, col=cols,lty=1,bty="n", cex=0.7)

	#matplot(R0test, t(start.difference.age), ylab=expression(hat(beta)[f]*" - "*hat(beta[a])), type="l", col=cols,lty=1, xlab=expression(R[0])); 
	#matplot(R0test, t(end.difference.age), ylab=expression(hat(beta)[f]*" - "*hat(beta[a])), type="l", col=cols,lty=1, xlab=expression(R[0])); 
	#matplot(R0test, t(max.difference.age), ylab=expression(hat(beta)[f]*" - "*hat(beta[a])), type="l", col=cols,lty=1, xlab=expression(R[0])); 
matplot(R0test, t(max.difference.self.age), ylab=expression(hat(beta)[a]*" - "*hat(beta[a1])), type="l", col=cols,lty=1, xlab=expression(R[0])); 

   par(mfrow=c(1,2), mar=c(4,5,2,1), bty="l")
   matplot(R0test,t(peak.incidence[2,,])/26,type="l", col=cols,lty=1, xlab=expression(R[0]), ylab="Time of peak/min incidence (yrs)")
   matplot(R0test,t(peak.incidence[3,,])/26,type="l", col=cols,lty=1, xlab=expression(R[0]), ylab="Time of peak/min incidence (yrs)")
legend("topright",legend= name.waifw, col=cols,lty=1,bty="n", cex=0.7)

   par(mfrow=c(2,2), mar=c(4,5,2,1), bty="l")
   matplot(R0test,t(incidence.at.start[2,,]),type="l", col=cols,lty=1, xlab=expression(R[0]), ylab="Incidence at start")
legend("topleft",legend= name.waifw, col=cols,lty=1,bty="n", cex=0.7)
   matplot(R0test,t(incidence.at.start[3,,]),type="l", col=cols,lty=1, xlab=expression(R[0]), ylab="Incidence at start")
   matplot(R0test,t(incidence.at.end[2,,]),type="l", col=cols,lty=1, xlab=expression(R[0]), ylab="Incidence at end")
   matplot(R0test,t(incidence.at.end[3,,]),type="l", col=cols,lty=1, xlab=expression(R[0]), ylab="Incidence at end")



   par(mfrow=c(1,2), mar=c(4,5,2,1), bty="l")
   matplot(R0test,t(organize.age.inc)/12,type="l", col=cols,lty=1, xlab=expression(R[0]), ylab="Age (start = line; end = dashed)", ylim=c(1,15))
   matplot(R0test,t(organize.age.inc.end)/12,type="l", col=cols,lty=3,add=TRUE)
   legend("topright",legend= name.waifw, col=cols,lty=1,bty="n", cex=0.7)
   matplot(R0test,t(organize.age.age)/12,type="l", col=cols,lty=1, xlab=expression(R[0]), ylab="Age (start = line; end = dashed)", ylim=c(1,15))
   matplot(R0test,t(organize.age.age.end)/12,type="l", col=cols,lty=3,add=TRUE)


   par(mfrow=c(1,2), mar=c(4,5,2,1), bty="l")
   matplot(R0test,t(organize.var.age.inc),type="l", col=cols,lty=1, xlab=expression(R[0]), ylab="Var age (start = line; end = dashed)")
   matplot(R0test,t(organize.var.age.inc.end),type="l", col=cols,lty=3,add=TRUE)
   legend("topright",legend= name.waifw, col=cols,lty=1,bty="n", cex=0.7)
   matplot(R0test,t(organize.var.age.age),type="l", col=cols,lty=1, xlab=expression(R[0]), ylab="Var age (start = line; end = dashed)")
   matplot(R0test,t(organize.var.age.age.end),type="l", col=cols,lty=3,add=TRUE)


par(mfrow=c(1,1))
matplot(t(incidences[incidences[,3]==14 & incidences[,1]==2,4:ncol(incidences)]), type="l", col=cols, lty=1, xlab="Time")



### PLOTS FOR THE PRESENTATION ############################################################################ 

### 1. Phase planes
par(mfrow=c(1,1), mar=c(5,5,3,3))
plot(t(susceptibles[incidences[,3]==14 & incidences[,1]==2 & incidences[,2]==1,4:ncol(incidences)]),t(incidences[incidences[,3]==14 & incidences[,1]==2 & incidences[,2]==1,4:ncol(incidences)]), 
	type="l", col=cols, lty=1, xlab="Number of susceptible individuals", ylab="Number of infected individuals", cex.lab=1.5, 
	xlim=range(c(susceptibles[incidences[,3]==14 & incidences[,1]==2,4:ncol(incidences)])), 
	ylim=range(c(incidences[incidences[,3]==14 & incidences[,1]==2,4:ncol(incidences)])),lwd=3)
abline(h=mean(incidences[incidences[,3]==14 & incidences[,1]==2,4]),lty=3, col="grey")

matplot(t(susceptibles[incidences[,3]==14 & incidences[,1]==2,4:ncol(incidences)]),t(incidences[incidences[,3]==14 & incidences[,1]==2,4:ncol(incidences)]), 
	type="l", col=cols, lty=1, xlab="Number of susceptible individuals", ylab="Number of infected individuals", cex.lab=1.5,lwd=3)
abline(h=mean(incidences[incidences[,3]==14 & incidences[,1]==2,4]),lty=3, col="grey")


### 2. Dissect dynamics ##########
chs.R0 <- 3 #corresonds to 14
chs.row.baseline <- c(1:length(index.waifw))[index.waifw == 1 & index.R0==R0test[chs.R0]]	
tx <- (1:length(chs.times))/26  #define the x axis
xlims <- c(9,20)

par(mfrow=c(3,1), mar=c(2,5,1,0.5),bty="l")
plot(tx ,as.numeric(incidences[chs.row.baseline[2], chs.times]), type="l", xlab="", ylab="# infected individuals", col=1,lwd=2,
		ylim=range(incidences[index.R0==R0test[chs.R0] & incidences[,1]==2,chs.times],na.rm=TRUE), cex.lab=1.5, xlim=xlims)
abline(v=tx[vacc.time],col="grey")  
for (k in 1:5) { 
chs.row <- c(1:length(index.waifw))[index.waifw == k & index.R0==R0test[chs.R0]]	
j <- 2; points(tx,as.numeric(incidences[chs.row[j], chs.times]),type="l",col=cols[k],lwd=2)
}
legend("topleft",legend= name.waifw, col=cols,lty=1,bty="n", cex=1,lwd=2)

plot(tx,as.numeric(store.beta.unity[chs.row.baseline[2], chs.times]), type="l", xlab="", ylab=expression(hat(beta)), col=1,lwd=2,
		ylim=range(store.beta.unity[index.R0==R0test[chs.R0],chs.times],na.rm=TRUE), cex.lab=1.5, xlim=xlims)
abline(v=tx[vacc.time],col="grey")  
for (k in 1:5) { 
chs.row <- c(1:length(index.waifw))[index.waifw == k & index.R0==R0test[chs.R0]]	
j <- 2; points(tx,as.numeric(store.beta.unity[chs.row[j], chs.times]),type="l",col=cols[k],lwd=2)
}


plot(tx,as.numeric(avg.age[chs.row.baseline[2], chs.times])/12, type="l", xlab="", ylab="Average age", col=1,lwd=2,
		ylim=range(avg.age[index.R0==R0test[chs.R0],chs.times],na.rm=TRUE)/12, cex.lab=1.5, xlim=xlims)
abline(v=tx[vacc.time],col="grey")  
abline(h=c(5,10,15),lty=3)
for (k in 1:5) { 
chs.row <- c(1:length(index.waifw))[index.waifw == k & index.R0==R0test[chs.R0]]	
j <- 2; points(tx,as.numeric(avg.age[chs.row[j], chs.times])/12,type="l",col=cols[k],lwd=2)
}




#### 3. One year later ##########
par(mfrow=c(1,2), mar=c(5,6,2,2))
matplot(R0test,t(max.difference.self.incidence[,]),type="b", col=cols,lty=1, xlab=expression(R[0]), 
		ylab=expression("Largest change in "*hat(beta)), cex.lab=1.5,lwd=3, pch=19)
matplot(R0test,t(difference.self.incidence.1yr[,]),type="b", col=cols,lty=1, xlab=expression(R[0]), 
		ylab=expression("Change in "*hat(beta)* " one year after vaccination"), cex.lab=1.5,lwd=3, pch=19)


### What happens to age (they move AWAY from the optimal and then the variance tightens up so moves in closer? But there is also how fast it is happening?)
layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))
matplot(R0test,t(difference.self.incidence.1yr[,]),type="b", col=cols,lty=1, xlab=expression(R[0]), 
		ylab=expression("Change in "*hat(beta)* " one year after vaccination"), cex.lab=1.5,lwd=3, pch=19)
matplot(R0test,t(organize.age.inc[,])/12,type="b", col=cols,lty=1, xlab=expression(R[0]), ylab=expression("Average age (start)"), cex.lab=1.5,lwd=3, pch=19, ylim=c(0,16))
legend("topright",legend= name.waifw, col=cols,lty=1,bty="n", cex=0.7,lwd=2)
matplot(R0test,t(organize.var.age.inc[,])/12,type="b", col=cols,lty=1, xlab=expression(R[0]), ylab=expression("Variance age (start)"), cex.lab=1.5,lwd=3, pch=19,log="y")


### What happens to age - show largest change and different at the end
layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))
matplot(R0test,t(max.difference.self.incidence[,]),type="b", col=cols,lty=1, xlab=expression(R[0]), 
		ylab=expression("Largest change in "*hat(beta)), cex.lab=1.5,lwd=3, pch=19)
matplot(R0test,t(organize.age.inc.end[,])/12,type="b", col=cols,lty=1, xlab=expression(R[0]), ylab=expression("Average age (end)"), cex.lab=1.5,lwd=3, pch=19, ylim=c(0,16))
legend("topright",legend= name.waifw, col=cols,lty=1,bty="n", cex=0.7,lwd=2)
matplot(R0test,t(organize.var.age.inc.end[,])/12,type="b", col=cols,lty=1, xlab=expression(R[0]), ylab=expression("Variance age (end)"), cex.lab=1.5,lwd=3, pch=19,log="y")



## change in age
matplot(R0test,t(organize.age.inc -organize.age.inc.end[,])/12,type="b", col=cols,lty=1, xlab=expression(R[0]), 
		ylab=expression("Average age (difference)"), cex.lab=1.5,lwd=3, pch=19)




### 3. Peak incidence
par(mfrow=c(1,2), mar=c(5,5,3,3), bty="l")
matplot(R0test,t(peak.incidence[2,,])/26,type="b", col=cols,lty=1, xlab=expression(R[0]), ylab="Time of peak/min incidence (yrs)", cex.lab=1.5,lwd=3, pch=19)
matplot(R0test,t(peak.incidence[3,,])/26,type="b", col=cols,lty=1, xlab=expression(R[0]), ylab="Time of peak/min incidence (yrs)", cex.lab=1.5,lwd=3, pch=19)
legend("topright",legend= name.waifw, col=cols,lty=1,bty="n", cex=0.7)



## no mega obvious pattern with the range 
plot(c(peak.incidence[2,,]),c(max.difference.self.incidence),col=cols, xlab="Time to peak", ylab=expression("Max change in "*hat(beta)), pch=19)
legend("topright",legend= name.waifw, col=cols,lty=1,bty="n", cex=0.7)
abline(h=0)

par(mfrow=c(1,1), mar=c(5,5,3,3), bty="l")
plot(c(peak.incidence[2,,]),c(difference.self.incidence.1yr[,]),col=cols, xlab="Time to peak", ylab=expression("Max change in "*hat(beta)), pch=19)
#legend("bottomright",legend= name.waifw, col=cols,lty=1,bty="n", cex=0.7)
abline(h=0, col="grey")

plot(c(cumsum.beta.1yr[,]),c(peak.incidence[2,,])/26,col=cols, ylab="Time to peak", xlab=expression("Cumulative sum "*hat(beta)*" up to 1 yr"), pch=19)
#legend("bottomright",legend= name.waifw, col=cols,lty=1,bty="n", cex=0.7)
abline(h=0, col="grey")

## can scale by R0 or /start.beta[,] or whatever you like but doesn't help in separating out the different waifws, always overlap




	
	### VACCINATE OR NOT... THE FIRST TIME-STEP SHOULD HAVE THE SAME UNITY BETA IF MATCHING WORKED THE SAME? 
	store.beta.unity1 <- read.csv(file=paste("/Users/cmetcalf/Documents/temp/unity.beta",TRUE,".csv", sep=""), header=T, row.names=NULL)
	store.beta.unity2 <- read.csv(file=paste("/Users/cmetcalf/Documents/temp/unity.beta",FALSE,".csv", sep=""), header=T, row.names=NULL)
	plot(store.beta.unity1[,400], store.beta.unity2[,400], col=cols[store.beta.unity1[,2]], pch=19)
	abline(0,1)


	### One intricate look at a focal one

        cols <- c(1,4,2) # blue is matched on incidence, red is matched on age; the first should be ignored - is just the uniform run for matching. 

	chs.R0 <- 3; k <- 2
	chs.row.baseline <- c(1:length(index.waifw))[index.waifw == 1 & index.R0==R0test[chs.R0]]	
 	chs.row <- c(1:length(index.waifw))[index.waifw == k & index.R0==R0test[chs.R0]]

	par(mfrow=c(1,3),bty="l",pty="m")

	## PHASE PLANE
   	plot(as.numeric(susceptibles[chs.row.baseline[2], chs.times]), as.numeric(incidences[chs.row.baseline[2], chs.times]), 
		type="l", xlab="susceptible", ylab="infected", col=1,lwd=2,
		ylim=range(incidences[chs.row,chs.times]), xlim=range(susceptibles[chs.row,chs.times]),main=name.waifw[index.waifw[chs.row[1]]]) 
        for (j in 2:length(chs.row)) points(as.numeric(susceptibles[chs.row[j], chs.times]), as.numeric(incidences[chs.row[j], chs.times]), 
		type="l",lwd=lwds[index.match[chs.row[j]]],col=cols[index.match[chs.row[j]]])
	## AGE
        plot(as.numeric(avg.age[chs.row.baseline[2], chs.times])/12, type="l", xlab="", ylab="Average Age", col=1,lwd=2,
		ylim=range(c(avg.age[chs.row.baseline[2], chs.times],avg.age[chs.row,chs.times]),na.rm=TRUE)/12)
        for (j in 2:length(chs.row)) points(as.numeric(avg.age[chs.row[j], chs.times])/12,type="l",lwd=lwds[index.match[chs.row[j]]],col=cols[index.match[chs.row[j]]])
	## UNITY 
        plot(as.numeric(store.beta.unity[chs.row.baseline[2], chs.times]), type="l", xlab="", ylab=expression(hat(beta)), col=1,lwd=2,
		#ylim=range(store.beta.unity[chs.row,chs.times]),main=name.waifw[index.waifw[chs.row[j]]])
		ylim=range(store.beta.unity[index.R0==R0test[chs.R0],chs.times],na.rm=TRUE))
        for (j in 2:length(chs.row)) points(as.numeric(store.beta.unity[chs.row[j], chs.times]),type="l",lwd=lwds[index.match[chs.row[j]]],col=cols[index.match[chs.row[j]]])
    

	## this is what the average age should be if Life expectancy is 50
	50/(R0test[chs.R0]-1)

	
	par(mfrow=c(1,1),bty="l",pty="m")	
	## check out to understand age beta_unity dip? 
	plot(as.numeric(incidences[chs.row.baseline[2], chs.times]), type="l",lwd=2, col="grey", ylim=c(0,470))
	points(as.numeric(incidences[chs.row[2], chs.times]), type="l",lwd=2, col="blue")
	points(as.numeric(incidences[chs.row[3], chs.times]), type="l",lwd=2, col="red")

	points(400*as.numeric(store.beta.unity[chs.row[2], chs.times])/max(as.numeric(store.beta.unity[chs.row[2], chs.times])), type="l",col="blue")
        points(400*as.numeric(store.beta.unity[chs.row[3], chs.times])/max(as.numeric(store.beta.unity[chs.row[3], chs.times])), type="l",col="red")

	points(400*as.numeric(avg.age[chs.row[2], chs.times])/max(as.numeric(avg.age[chs.row[2], chs.times])), type="l",col="blue", lty=3)
	points(400*as.numeric(avg.age[chs.row[3], chs.times])/max(as.numeric(avg.age[chs.row[3], chs.times])), type="l",col="red", lty=3)
	### Flat incidence, flat average age.... ???
	



chs.R0 <- 4
chs.row <- c(1:length(index.waifw))[index.waifw == 4 & index.R0==R0test[chs.R0]]	
chs.row.baseline <- c(1:length(index.waifw))[index.waifw == 1 & index.R0==R0test[chs.R0]]	


par(mfrow=c(4,1),mar=c(2,4,1,2))
plot(as.numeric(store.beta.unity[chs.row[j], chs.times]),  as.numeric(incidences[chs.row[j], chs.times]),type="l", xlab=expression(beta),ylab="Infected")
plot(as.numeric(store.beta.unity[chs.row[j], chs.times]),  as.numeric(susceptibles[chs.row[j], chs.times]),type="l", xlab=expression(beta),ylab="Susceptible")
plot(as.numeric(store.beta.unity[chs.row[j], chs.times]),  1/(as.numeric(incidences[chs.row[j], chs.times]^0.97)*as.numeric(susceptibles[chs.row[j], chs.times])),
	type="l", xlab=expression(beta),ylab="1/((Infected^gamma)*Susceptible)")
plot(as.numeric(store.beta.unity[chs.row[j], chs.times]),  as.numeric(avg.age[chs.row[j], chs.times]),type="l", xlab=expression(beta),ylab="Avg age")

par(mfrow=c(5,1),mar=c(2,4,1,2))
maxB <- which(as.numeric(store.beta.unity[chs.row[j], chs.times])==max(as.numeric(store.beta.unity[chs.row[j], chs.times])))
plot(as.numeric(store.beta.unity[chs.row[j], chs.times]),type="l", xlab=expression(time),ylab="beta", 
	ylim=range(as.numeric(store.beta.unity[chs.row[j], chs.times]),as.numeric(store.beta.unity[chs.row.baseline[j], chs.times]))); abline(v=maxB)
points(as.numeric(store.beta.unity[chs.row.baseline[j], chs.times]),type="l", col="red")
plot(as.numeric(susceptibles[chs.row[j], chs.times]),type="l", xlab=expression(time),ylab="susceptible"); abline(v=maxB)
points(as.numeric(susceptibles[chs.row.baseline[j], chs.times]),type="l", col="red")
plot(as.numeric(incidences[chs.row[j], chs.times]),type="l", xlab=expression(time),ylab="Infected"); abline(v=maxB)
points(as.numeric(incidences[chs.row.baseline[j], chs.times]),type="l", col="red")
plot(1/(as.numeric(incidences[chs.row[j], chs.times]^0.97)*as.numeric(susceptibles[chs.row[j], chs.times])),type="l", ylab=expression(1/(SIalph)),xlab=""); abline(v=maxB)
plot(as.numeric(avg.age[chs.row[j], chs.times]),type="l", xlab=expression(time),ylab="Age"); abline(v=maxB)
points(as.numeric(avg.age[chs.row.baseline[j], chs.times]),type="l", col="red")


}


doAfew <- function(){
  
  do.plot <- FALSE
  remove.vacc <- FALSE #TRUE# #
  do.unity <- TRUE
  age.classes <- c(1:120,seq(180,600,by=60))
 
  background <- 0.001

  ## series from flat to peaky to diagonal to polymod
  waifw1 <- matrix(1,length(c(1:120,seq(180,600,by=60))/12),length(c(1:120,seq(180,600,by=60))/12))
  waifw2 <- get.smooth.WAIFW(age.class.boundries=c(1:120,seq(180,600,by=60))/12, mu=5,sig=0.2, gam=0.05, delta= background)
  waifw3 <- get.smooth.WAIFW(age.class.boundries=c(1:120,seq(180,600,by=60))/12, mu=10,sig=0.1, gam=0.01, delta= background)
  waifw4 <- get.smooth.WAIFW(age.class.boundries=c(1:120,seq(180,600,by=60))/12, mu=7,sig=0.6, gam=0.05, delta= background) 
  waifw5 <- get.polymod.WAIFW(age.class.boundries=c(1:120,seq(180,600,by=60))/12, do.touch=TRUE)
 

  #waifw2 <- get.smooth.WAIFW(age.class.boundries=c(1:120,seq(180,600,by=60))/12, mu=6,sig=0.5, gam=0.05, delta=0.1) 
  #waifw3 <- get.smooth.WAIFW(age.class.boundries=c(1:120,seq(180,600,by=60))/12, mu=30,sig=1, gam=0.08, delta=0.1) 

  
  z10 <- aFullComparisonSimple(R0=10, remove.vacc= remove.vacc, waifw.test=waifw1, full.plot=do.plot,index=1, do.unity= do.unity)
  z15 <- aFullComparisonSimple(R0=12, remove.vacc= remove.vacc, waifw.test=waifw1, full.plot=do.plot,index=1, do.unity= do.unity)
  z20 <- aFullComparisonSimple(R0=14, remove.vacc= remove.vacc, waifw.test=waifw1, full.plot=do.plot,index=1, do.unity= do.unity)
  z25 <- aFullComparisonSimple(R0=16, remove.vacc= remove.vacc, waifw.test=waifw1, full.plot=do.plot,index=1, do.unity= do.unity)
  z30 <- aFullComparisonSimple(R0=18, remove.vacc= remove.vacc, waifw.test=waifw1, full.plot=do.plot,index=1, do.unity= do.unity)
  z35 <- aFullComparisonSimple(R0=20, remove.vacc= remove.vacc, waifw.test=waifw1, full.plot=do.plot,index=1, do.unity= do.unity)
  z40 <- aFullComparisonSimple(R0=22, remove.vacc= remove.vacc, waifw.test=waifw1, full.plot=do.plot,index=1, do.unity= do.unity)
  z45 <- aFullComparisonSimple(R0=24, remove.vacc= remove.vacc, waifw.test=waifw1, full.plot=do.plot,index=1, do.unity= do.unity)
 
  z10a <- aFullComparisonSimple(R0=10, remove.vacc= remove.vacc, waifw.test=waifw2, full.plot=do.plot,index=2, do.unity= do.unity)
  z15a <- aFullComparisonSimple(R0=12, remove.vacc= remove.vacc, waifw.test=waifw2, full.plot=do.plot,index=2, do.unity= do.unity)
  z20a <- aFullComparisonSimple(R0=14, remove.vacc= remove.vacc, waifw.test=waifw2, full.plot=do.plot,index=2, do.unity= do.unity)
  z25a <- aFullComparisonSimple(R0=16, remove.vacc= remove.vacc, waifw.test=waifw2, full.plot=do.plot,index=2, do.unity= do.unity)
  z30a <- aFullComparisonSimple(R0=18, remove.vacc= remove.vacc, waifw.test=waifw2, full.plot=do.plot,index=2, do.unity= do.unity)
  z35a <- aFullComparisonSimple(R0=20, remove.vacc= remove.vacc, waifw.test=waifw2, full.plot=do.plot,index=2, do.unity= do.unity)
  z40a <- aFullComparisonSimple(R0=22, remove.vacc= remove.vacc, waifw.test=waifw2, full.plot=do.plot,index=2, do.unity= do.unity)
  z45a <- aFullComparisonSimple(R0=24, remove.vacc= remove.vacc, waifw.test=waifw2, full.plot=do.plot,index=2, do.unity= do.unity)

  z10b <- aFullComparisonSimple(R0=10, remove.vacc= remove.vacc, waifw.test=waifw3, full.plot=do.plot,index=3, do.unity= do.unity)
  z15b <- aFullComparisonSimple(R0=12, remove.vacc= remove.vacc, waifw.test=waifw3, full.plot=do.plot,index=3, do.unity= do.unity)
  z20b <- aFullComparisonSimple(R0=14, remove.vacc= remove.vacc, waifw.test=waifw3, full.plot=do.plot,index=3, do.unity= do.unity)
  z25b <- aFullComparisonSimple(R0=16, remove.vacc= remove.vacc, waifw.test=waifw3, full.plot=do.plot,index=3, do.unity= do.unity)
  z30b <- aFullComparisonSimple(R0=18, remove.vacc= remove.vacc, waifw.test=waifw3, full.plot=do.plot,index=3, do.unity= do.unity)
  z35b <- aFullComparisonSimple(R0=20, remove.vacc= remove.vacc, waifw.test=waifw3, full.plot=do.plot,index=3, do.unity= do.unity)
  z40b <- aFullComparisonSimple(R0=22, remove.vacc= remove.vacc, waifw.test=waifw3, full.plot=do.plot,index=3, do.unity= do.unity)
  z45b <- aFullComparisonSimple(R0=24, remove.vacc= remove.vacc, waifw.test=waifw3, full.plot=do.plot,index=3, do.unity= do.unity)
 
  z10c <- aFullComparisonSimple(R0=10, remove.vacc= remove.vacc, waifw.test=waifw4, full.plot=do.plot,index=4, do.unity= do.unity)
  z15c <- aFullComparisonSimple(R0=12, remove.vacc= remove.vacc, waifw.test=waifw4, full.plot=do.plot,index=4, do.unity= do.unity)
  z20c <- aFullComparisonSimple(R0=14, remove.vacc= remove.vacc, waifw.test=waifw4, full.plot=do.plot,index=4, do.unity= do.unity)
  z25c <- aFullComparisonSimple(R0=16, remove.vacc= remove.vacc, waifw.test=waifw4, full.plot=do.plot,index=4, do.unity= do.unity)
  z30c <- aFullComparisonSimple(R0=18, remove.vacc= remove.vacc, waifw.test=waifw4, full.plot=do.plot,index=4, do.unity= do.unity)
  z35c <- aFullComparisonSimple(R0=20, remove.vacc= remove.vacc, waifw.test=waifw4, full.plot=do.plot,index=4, do.unity= do.unity)
  z40c <- aFullComparisonSimple(R0=22, remove.vacc= remove.vacc, waifw.test=waifw4, full.plot=do.plot,index=4, do.unity= do.unity)
  z45c <- aFullComparisonSimple(R0=24, remove.vacc= remove.vacc, waifw.test=waifw4, full.plot=do.plot,index=4, do.unity= do.unity)
  
  z10d <- aFullComparisonSimple(R0=10, remove.vacc= remove.vacc, waifw.test=waifw5, full.plot=do.plot,index=5, do.unity= do.unity)
  z15d <- aFullComparisonSimple(R0=12, remove.vacc= remove.vacc, waifw.test=waifw5, full.plot=do.plot,index=5, do.unity= do.unity)
  z20d <- aFullComparisonSimple(R0=14, remove.vacc= remove.vacc, waifw.test=waifw5, full.plot=do.plot,index=5, do.unity= do.unity)
  z25d <- aFullComparisonSimple(R0=16, remove.vacc= remove.vacc, waifw.test=waifw5, full.plot=do.plot,index=5, do.unity= do.unity)
  z30d <- aFullComparisonSimple(R0=18, remove.vacc= remove.vacc, waifw.test=waifw5, full.plot=do.plot,index=5, do.unity= do.unity)
  z35d <- aFullComparisonSimple(R0=20, remove.vacc= remove.vacc, waifw.test=waifw5, full.plot=do.plot,index=5, do.unity= do.unity)
  z40d <- aFullComparisonSimple(R0=22, remove.vacc= remove.vacc, waifw.test=waifw5, full.plot=do.plot,index=5, do.unity= do.unity)
  z45d <- aFullComparisonSimple(R0=24, remove.vacc= remove.vacc, waifw.test=waifw5, full.plot=do.plot,index=5, do.unity= do.unity)
  
  
  allZ <- rbind(z10,z15,z20,z25,z30,z35,z40,z45)
  allZa <- rbind(z10a,z15a,z20a,z25a,z30a,z35a,z40a,z45a)
  allZb <- rbind(z10b,z15b,z20b,z25b,z30b,z35b,z40b,z45b)
  allZc <- rbind(z10c,z15c,z20c,z25c,z30c,z35c,z40c,z45c)
  allZd <- rbind(z10d,z15c,z20d,z25d,z30d,z35d,z40d,z45d)
  
  cols <- rep(colorRampPalette(c(rgb(0,0,1,0.8), rgb(1,0,0,0.8)), alpha = TRUE)(8), each =3)
  pchs <- rep(c(19,15,3),8)




  ### 1. Dissecting the scaling ######################################
  ## Plot out the waifws, and below, plot the relationship between the target R0 and max waifw.  
 
  par(mfrow=c(1,5),pty="s")
  image(age.classes/12,age.classes/12, waifw1, xlab="age", ylab="age", main="",col = hcl.colors(15, alpha=0.4, rev = TRUE),zlim=c(waifw1[1,1],1000)); abline(h=5, lty=3); abline(v=5,lty=3)
  image(age.classes/12,age.classes/12, waifw2, xlab="age", ylab="age", main="",col = hcl.colors(15, alpha=0.4, rev = TRUE)); abline(h=5, lty=3); abline(v=5,lty=3)
  image(age.classes/12,age.classes/12, waifw3, xlab="age", ylab="age", main="",col = hcl.colors(15, alpha=0.4, rev = TRUE)); abline(h=5, lty=3); abline(v=5,lty=3)
  image(age.classes/12,age.classes/12, waifw4, xlab="age", ylab="age", main="",col = hcl.colors(15, alpha=0.4, rev = TRUE)); abline(h=5, lty=3); abline(v=5,lty=3)
  image(age.classes/12,age.classes/12, waifw5, xlab="age", ylab="age", main="Polymod",col = hcl.colors(15, alpha=0.4, rev = TRUE)); abline(h=5, lty=3); abline(v=5,lty=3)

  for (j in 1:5) { 
    if (j==1) testX <- allZ;  if (j==2) testX <- allZa;  if (j==3) testX <- allZb;  if (j==4) testX <- allZc;  if (j==5) testX <- allZd; 
    plot(testX$target.R0[testX$type.match=="incidence"],testX$max.waifw[testX$type.match=="incidence"], type="b", ylim=range(c(allZ$max.waifw,allZa$max.waifw,allZb$max.waifw,allZc$max.waifw)),
		pch=19,col="lightgrey", xlab=expression("Target "*R[0]), ylab="Max transmission")
    points(testX$target.R0[testX$type.match=="age"],testX$max.waifw[testX$type.match=="age"], type="b",pch=19,col="darkgrey")
    points(testX$target.R0[testX$type.match=="none"],testX$max.waifw[testX$type.match=="none"], type="b",pch=19, col="black")
  }
  legend("topleft",legend=c("baseline","matched to incidence", "matched to age"), lty=1,col=c("black","lightgrey","darkgrey"),bty="n")

  ## sanity check - max waifw for the flat should have a linear relationship with the target R0 with slope of 1 scale by the population size
  #fit <- lm(allZb$max.waifw[allZ$type.match=="none"]~allZb$target.R0[allZ$type.match=="none"])
  ## 500000*fit$coeff[2]
  ## to match on incidence, you seem to have to go higher (slope is steeper); to match on age, slope is even steeper
  ## General pattern is one of simple increase 




  ### 2. Exploring impacts on dynamics ######################################
  ## Note that as R0 goes up, incidence pre-vacc goes up; and average age goes down. 
  ## these projections work for vaccine release

  ## effects on timing and height of peak? - columns are different WAIFWs, pch is i) flat for that R0, ii) matched on incidnece, ii) then age (pch=19,15,3); colours are R0
  par(mfrow=c(2,5))
  plot(allZ$target.R0,allZ$tpeak/26,xlab="Target R0", ylab="Time to peak (years)", col=cols, pch=pchs, ylim=range(c(allZ$tpeak,allZa$tpeak,allZb$tpeak,allZc$tpeak,allZd$tpeak)/26))
  plot(allZa$target.R0,allZa$tpeak/26,xlab="Target R0", ylab="Time to peak (years)", col=cols, pch=pchs, ylim=range(c(allZ$tpeak,allZa$tpeak,allZb$tpeak,allZc$tpeak,allZd$tpeak)/26))
  plot(allZa$target.R0,allZb$tpeak/26,xlab="Target R0", ylab="Time to peak (years)", col=cols, pch=pchs, ylim=range(c(allZ$tpeak,allZa$tpeak,allZb$tpeak,allZc$tpeak,allZd$tpeak)/26))
  plot(allZa$target.R0,allZc$tpeak/26,xlab="Target R0", ylab="Time to peak (years)", col=cols, pch=pchs, ylim=range(c(allZ$tpeak,allZa$tpeak,allZb$tpeak,allZc$tpeak,allZd$tpeak)/26))
  plot(allZa$target.R0,allZd$tpeak/26,xlab="Target R0", ylab="Time to peak (years)", col=cols, pch=pchs, ylim=range(c(allZ$tpeak,allZa$tpeak,allZb$tpeak,allZc$tpeak,allZd$tpeak)/26))



  #plot(allZ$target.R0,allZ$incidence.peak, log="y",xlab="Target R0", ylab="Max peak", col=cols, pch=pchs)
  #plot(allZa$target.R0,allZa$incidence.peak, log="y",xlab="Target R0", ylab="Max peak", col=cols, pch=pchs)
  #plot(allZa$target.R0,allZb$incidence.peak, log="y",xlab="Target R0", ylab="Max peak", col=cols, pch=pchs)
  #plot(allZa$target.R0,allZc$incidence.peak, log="y",xlab="Target R0", ylab="Max peak", col=cols, pch=pchs)
  #plot(allZa$target.R0,allZd$incidence.peak, log="y",xlab="Target R0", ylab="Max peak", col=cols, pch=pchs)


  ## Vaccine release
  ## 1st column; as R0 increases delay goes down, and max peak also goes cos less time susceptible accumulatiion (?paradox of enrichment) [very low R0 delay and peak goes down?] 
  ## 2nd column: delay first increases (less susceptibles of the right age) then decreases (susceptibles depleted generally); max peak goes down first (why? taking you longer but max peak transmission scaling offsets this?)
  ## 3rd column: we're back to down and down
  ## 4th column: delay goes down but peak goes up, rather like the 2nd half of second.  


  ylims <- c(4,40)
 
  ## same as second, but rescale by incidence, multiply by 26 because the incidence pre-peak was annual
  plot(allZ$target.R0,26*allZ$incidence.peak/allZ$incidence.pre.event, log="y",
       xlab="Target R0", ylab="Relative peak", col=cols, pch=pchs, ylim=ylims)
  plot(allZa$target.R0,26*allZa$incidence.peak/allZa$incidence.pre.event, log="y",
       xlab="Target R0", ylab="Relative peak", col=cols, pch=pchs, ylim=ylims)
  plot(allZb$target.R0,26*allZb$incidence.peak/allZb$incidence.pre.event, log="y",
       xlab="Target R0", ylab="Relative peak", col=cols, pch=pchs, ylim=ylims)
  plot(allZc$target.R0,26*allZc$incidence.peak/allZc$incidence.pre.event, log="y",
       xlab="Target R0", ylab="Relative peak", col=cols, pch=pchs, ylim=ylims)
  plot(allZd$target.R0,26*allZd$incidence.peak/allZc$incidence.pre.event, log="y",
       xlab="Target R0", ylab="Relative peak", col=cols, pch=pchs, ylim=ylims)



  ylims <- c(1,12)

  fit.uniform.s <- lm(I(allZ$age.after.five.years/12)~I(allZ$age.pre.event/12),subset=(allZ$type.waifw=="flat"))
  fit.uniform.l <- lm(I(allZ$age.at.end/12)~I(allZ$age.pre.event/12),subset=(allZ$type.waifw=="flat"))

  ##
  par(mfrow=c(1,5))
  plot(allZ$age.pre.event/12 ,allZ$age.after.five.years/12, col=cols, pch=pchs,xlab="Age before", ylab="Age after", ylim=ylims)
  points(allZ$age.pre.event/12 ,allZ$age.at.end/12, col=cols, pch=c(1,0,4))

  legend("topleft",legend=c("5 yr flat","5 yr incidence match","5 yr age match","long flat","long incidence match","long age match"), pch=c(pchs[1:3],c(1,0,4)),bty="n")

  abline(0,1,lty=1,lwd=0.5); abline(fit.uniform.s,lty=3); abline(fit.uniform.l,lty=3)
  plot(allZa$age.pre.event/12 ,allZa$age.after.five.years/12, col=cols, pch=pchs,xlab="Age before", ylab="Age after", ylim=ylims)
  points(allZa$age.pre.event/12 ,allZa$age.at.end/12, col=cols, pch=c(1,0,4))
  abline(0,1,lty=1,lwd=0.5)
abline(fit.uniform.s,lty=3); abline(fit.uniform.l,lty=3)
  plot(allZb$age.pre.event/12,allZb$age.after.five.years/12, col=cols, pch=pchs,xlab="Age before", ylab="Age after", ylim=ylims)
  points(allZb$age.pre.event/12,allZb$age.at.end/12, col=cols, pch=c(1,0,4))
  abline(0,1,lty=1,lwd=0.5)
abline(fit.uniform.s,lty=3); abline(fit.uniform.l,lty=3)
  plot(allZc$age.pre.event/12,allZc$age.after.five.years/12, col=cols, pch=pchs,xlab="Age before", ylab="Age after", ylim=ylims)
  points(allZc$age.pre.event/12,allZc$age.at.end/12, col=cols, pch=c(1,0,4))
  abline(0,1,lty=1,lwd=0.5)
abline(fit.uniform.s,lty=3); abline(fit.uniform.l,lty=3)
  plot(allZd$age.pre.event/12,allZd$age.after.five.years/12, col=cols, pch=pchs,xlab="Age before", ylab="Age after", ylim=ylims)
  points(allZd$age.pre.event/12,allZd$age.at.end/12, col=cols, pch=c(1,0,4))
  abline(0,1,lty=1,lwd=0.5)
abline(fit.uniform.s,lty=3); abline(fit.uniform.l,lty=3)



## the lines corresponding to 'fit.uniform' tells you were would go in short and longer term. 
## the super spikey ones stay about the same, so between the uniform lines, diagonal looks like uniform, polymod a bit between the two. 

## is it at all possible to predict from first principles? 



  ### 3. Checks on how it is working ######################################



  #how good is the matching?  
  idx <- seq(1,24,by=3)
  par(mfrow=c(2,2))
   j<-1 #compare matched on incidence
   plot(allZ$incidence.pre.event[idx], allZ$incidence.pre.event[idx+j],pch=19, xlab="Incidence in flat", ylab="Incidence in model", 
		ylim=range(c(allZ$incidence.pre.event[idx+j],allZa$incidence.pre.event[idx+j],allZb$incidence.pre.event[idx+j],allZc$incidence.pre.event[idx+j]),na.rm=TRUE)); abline(0,1)
   points(allZa$incidence.pre.event[idx], allZa$incidence.pre.event[idx+j],pch=19,col="purple")
   points(allZb$incidence.pre.event[idx], allZb$incidence.pre.event[idx+j],pch=19,col="orange")
   points(allZc$incidence.pre.event[idx], allZc$incidence.pre.event[idx+j],pch=19,col="green")
   points(allZd$incidence.pre.event[idx], allZd$incidence.pre.event[idx+j],pch=19,col="pink")

   j<-2 #compare matched on age
   plot(allZ$age.pre.event[idx]/12, allZ$age.pre.event[idx+j]/12, xlab="Age in flat", ylab="Age in model",pch=19, ylim=c(1,15)); abline(0,1)
   points(allZa$age.pre.event[idx]/12, allZa$age.pre.event[idx+j]/12, pch=19,col="purple"); 
   points(allZb$age.pre.event[idx]/12, allZb$age.pre.event[idx+j]/12, pch=19,col="orange"); 
   points(allZc$age.pre.event[idx]/12, allZc$age.pre.event[idx+j]/12, pch=19,col="green"); 
   points(allZc$age.pre.event[idx]/12, allZc$age.pre.event[idx+j]/12, pch=19,col="pink"); 
 
  


  ## effects on age at the end? 
  par(mfrow=c(2,4))
  plot(allZ$target.R0,allZ$age.at.end/12, log="y",
       xlab="Target R0", ylab="Age at end", col=cols, pch=pchs)
  plot(allZ$target.R0,allZa$age.at.end/12, log="y",
       xlab="Target R0", ylab="Age at end", col=cols, pch=pchs)
  plot(allZ$target.R0,allZb$age.at.end/12, log="y",
       xlab="Target R0", ylab="Age at end", col=cols, pch=pchs)
  plot(allZ$target.R0,allZc$age.at.end/12, log="y",
       xlab="Target R0", ylab="Age at end", col=cols, pch=pchs)
  
  plot(allZ$target.R0,allZ$age.pre.event/12, log="y",
       xlab="Target R0", ylab="Age at start", col=cols, pch=pchs)
  plot(allZ$target.R0,allZa$age.pre.event/12, log="y",
       xlab="Target R0", ylab="Age at start", col=cols, pch=pchs)
  plot(allZ$target.R0,allZb$age.pre.event/12, log="y",
       xlab="Target R0", ylab="Age at start", col=cols, pch=pchs)
  plot(allZ$target.R0,allZc$age.pre.event/12, log="y",
       xlab="Target R0", ylab="Age at start", col=cols, pch=pchs)
  
  par(mfrow=c(2,2))
  plot(allZ$age.pre.event/12,allZ$age.at.end/12, #log="y",
       xlab="Age at start", ylab="Age after five years", col=cols, pch=pchs); abline(0,1)
  plot(allZa$age.pre.event/12,allZa$age.at.end/12, #log="y",
       xlab="Age at start", ylab="Age after five years", col=cols, pch=pchs); abline(0,1)
  plot(allZb$age.pre.event/12,allZb$age.at.end/12, #log="y",
       xlab="Age at start", ylab="Age after five years", col=cols, pch=pchs); abline(0,1)
  plot(allZc$age.pre.event/12,allZc$age.at.end/12, #log="y",
       xlab="Age at start", ylab="Age after five years", col=cols, pch=pchs); abline(0,1)
  

  par(mfrow=c(2,2))
  plot(allZ$age.pre.event/12,allZ$age.after.five.years/12, #log="y",
       xlab="Age at start", ylab="Age after 5 years", col=cols, pch=pchs); abline(0,1); abline(1,1)
  plot(allZa$age.pre.event/12,allZa$age.after.five.years/12, #log="y",
       xlab="Age at start", ylab="Age after 5 years", col=cols, pch=pchs); abline(0,1); abline(1,1)
  plot(allZb$age.pre.event/12,allZb$age.after.five.years/12, #log="y",
       xlab="Age at start", ylab="Age after 5 years", col=cols, pch=pchs); abline(0,1); abline(1,1)
  plot(allZc$age.pre.event/12,allZc$age.after.five.years/12, #log="y",
       xlab="Age at start", ylab="Age after 5 years", col=cols, pch=pchs); abline(0,1); abline(1,1)
  
    
    #allZ$age.after.five.years
  
  
  # look at rel with age vs. predicted effect of R0 to make
  par(mfrow=c(2,4))
  plot(50/(allZ$target.R0),allZ$age.pre.event/12,pch=pchs); abline(0,1)
  plot(50/(allZa$target.R0),allZa$age.pre.event/12,pch=pchs); abline(0,1)
  plot(50/(allZb$target.R0),allZb$age.pre.event/12,pch=pchs); abline(0,1)
  plot(50/(allZb$target.R0),allZc$age.pre.event/12,pch=pchs); abline(0,1)
  plot(50/(allZ$target.R0),allZ$age.at.end/12,pch=pchs); abline(0,1)
  plot(50/(allZa$target.R0),allZa$age.at.end/12,pch=pchs); abline(0,1)
  plot(50/(allZb$target.R0),allZb$age.at.end/12,pch=pchs); abline(0,1)
  plot(50/(allZb$target.R0),allZc$age.at.end/12,pch=pchs); abline(0,1)
  
  ## is the size of the suceptible pool unchanging? (Fine and Clarkson conjecture )
  par(mfrow=c(1,4))
  plot(50/(allZ$target.R0-1),allZ$susceptible.pre.event,pch=pchs, ylim=range(c(allZ$susceptible.pre.event,allZ$susceptible.at.end)),log="y", ylab="Susceptibles"); abline(0,1)
  points(50/(allZ$target.R0-1),allZ$susceptible.at.end,pch=pchs, col=2);
  plot(50/(allZa$target.R0-1),allZa$susceptible.pre.event,pch=pchs, ylim=range(c(allZa$susceptible.pre.event,allZa$susceptible.at.end)),log="y", ylab="Susceptibles"); abline(0,1)
  points(50/(allZa$target.R0-1),allZa$susceptible.at.end,pch=pchs, col=2); 
  plot(50/(allZb$target.R0-1),allZb$susceptible.pre.event,pch=pchs, ylim=range(c(allZb$susceptible.pre.event,allZb$susceptible.at.end)),log="y", ylab="Susceptibles"); abline(0,1)
  points(50/(allZb$target.R0-1),allZb$susceptible.at.end,pch=pchs, col=2); 
  plot(50/(allZb$target.R0-1),allZc$susceptible.pre.event,pch=pchs, ylim=range(c(allZc$susceptible.pre.event,allZc$susceptible.at.end)),log="y", ylab="Susceptibles"); abline(0,1)
  points(50/(allZb$target.R0-1),allZc$susceptible.at.end,pch=pchs, col=2); 
  
  ## this set is a good spectrum!
  par(mfrow=c(4,1))
  plot(allZ$R0~allZ$target.R0,pch=pchs); abline(0,1)
  plot(allZa$R0~allZa$target.R0,pch=pchs); abline(0,1)
  plot(allZb$R0~allZb$target.R0,pch=pchs); abline(0,1)
  plot(allZc$R0~allZc$target.R0,pch=pchs); abline(0,1)
  
par(mfrow=c(2,2))  
plot(allZ$R0,(allZ$age.pre.event-allZ$age.after.five.years)/12,pch=pchs)
plot(allZa$R0,(allZa$age.pre.event-allZa$age.after.five.years)/12,pch=pchs)
plot(allZb$R0,(allZb$age.pre.event-allZb$age.after.five.years)/12,pch=pchs)
plot(allZc$R0,(allZc$age.pre.event-allZc$age.after.five.years)/12,pch=pchs)
##only the super peaked gives you change of <1 year in age for the super young

par(mfrow=c(2,2))
plot((allZ$age.pre.event-allZ$age.after.five.years)[idx]/12,(allZ$age.pre.event-allZ$age.after.five.years)[idx+1]/12,pch=pchs[idx+1], col=cols[idx+1], ylim=c(-3,0), xlim=c(-3,0)); abline(0,1); abline(h=c(0:-2),lty=3); abline(v=c(0:-2),lty=3)
points((allZ$age.pre.event-allZ$age.after.five.years)[idx]/12,(allZ$age.pre.event-allZ$age.after.five.years)[idx+2]/12,pch=pchs[idx+2], col=cols[idx+12])

plot((allZa$age.pre.event-allZa$age.after.five.years)[idx]/12,(allZa$age.pre.event-allZa$age.after.five.years)[idx+1]/12,pch=pchs[idx+1], col=cols[idx+1], ylim=c(-3,0), xlim=c(-3,0)); abline(0,1); abline(h=c(0:-2),lty=3); abline(v=c(0:-2),lty=3)
points((allZa$age.pre.event-allZa$age.after.five.years)[idx]/12,(allZa$age.pre.event-allZa$age.after.five.years)[idx+2]/12,pch=pchs[idx+2], col=cols[idx+12])

plot((allZb$age.pre.event-allZb$age.after.five.years)[idx]/12,(allZb$age.pre.event-allZb$age.after.five.years)[idx+1]/12,pch=pchs[idx+1], col=cols[idx+1], ylim=c(-3,0), xlim=c(-3,0)); abline(0,1); abline(h=c(0:-2),lty=3); abline(v=c(0:-2),lty=3)
points((allZb$age.pre.event-allZb$age.after.five.years)[idx]/12,(allZb$age.pre.event-allZb$age.after.five.years)[idx+2]/12,pch=pchs[idx+2], col=cols[idx+12])

plot((allZc$age.pre.event-allZc$age.after.five.years)[idx]/12,(allZc$age.pre.event-allZc$age.after.five.years)[idx+1]/12,pch=pchs[idx+1], col=cols[idx+1], ylim=c(-3,0), xlim=c(-3,0)); abline(0,1); abline(h=c(0:-2),lty=3); abline(v=c(0:-2),lty=3)
points((allZc$age.pre.event-allZc$age.after.five.years)[idx]/12,(allZc$age.pre.event-allZc$age.after.five.years)[idx+2]/12,pch=pchs[idx+2], col=cols[idx+12])
## want to be above the zero one line (smaller age gap then uniform )
## suggests that diagonal (middle) not great
  

###########
#### Next thing is to pick out an intermediate level of R0, repeat with seasonality of 0.2 and show what that does, and then bring in Unity? 


}

##https://bmcpublichealth.biomedcentral.com/articles/10.1186/1471-2458-8-338
##Coverage gradually improved from approximately 50% during the 1970s to 86% when MMR vaccine replaced single antigen vaccine in 1988, and reached 92% in 1995.

## Age / R0 vaccination is tremendously perturbed if populations growing / shrinking (as expected)
## This might make 

## Birth rate goes from 18.5 to 11.1 per 1000



### TODO - check that not matching on the year AFTER the vaccination change? ##