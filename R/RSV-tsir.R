
simulateSIR <- function(R0=10,alpha=0.2,alpha.births=0,offset.births=0,startpop=100000,birth.rate=0.05,intro.rate=0.01,nsim=10,burn.in=24*115,Tmax=24*120, stochastic=FALSE) {  
 
    require(MASS)

    # seasonality and R0 (i.e., seasonal flucutations in transmission)
    season.force <- R0*(1+alpha*cos(2*pi*(1:52)/52))
                    
    # season.index
    season.index <- rep(1:52,Tmax/52)

    # get births per time-step (assume the birth rate provided is annual)
    births <- startpop*birth.rate*(1+alpha.births*cos(2*pi*(1:52+offset.births)/52))	

    # storage
    storeI <- storeS <- matrix(NA,nsim,Tmax)
    storeS[,1] <- startpop 
    storeI[,1] <- 10

    # loop
    for (j in 2:Tmax) {
        lambda <- season.force[season.index[j]]*(storeI[,j-1]^0.97)*storeS[,j-1]/startpop   		#calculate FOI for this season.
        if (stochastic) storeI[,j] <-rnegbin(nsim, lambda, pmax(storeI[,j - 1], 1))+         	        #new infections
            			rbinom(nsim,1,intro.rate)   			 	        #infected immigrants from elsewhere

	if (!stochastic) storeI[,j] <- lambda

        if (j<burn.in) storeI[,j-1] <- pmax(storeI[,j-1],10)  			       		#during burnin, make sure no extinction
        storeS[,j] <- pmax(storeS[,j-1]-storeI[,j]+births[season.index[j]],1) 	       		#make sure you have at least one susceptible
    }

	return(list(storeI=storeI,storeS=storeS,burn.in=burn.in,Tmax=Tmax))
}




## trial run without any seasonality in births
a1 <- simulateSIR(R0=7,alpha=0.3,alpha.births=0,offset.births=0,startpop=22709527,birth.rate= 0.0003212749,intro.rate=0.005,nsim=5,burn.in=52*115,Tmax=52*120)

## trial run with seasonality in births that also has an offset, so that they are different
a2 <- simulateSIR(R0=7,alpha=0.3,alpha.births=10,offset.births=26,startpop=22709527,birth.rate= 0.0003212749,intro.rate=0.005,nsim=5,burn.in=52*115,Tmax=52*120)

## plot them both together 
cols <- colorRampPalette(c("blue", "red"))(5) 
par(mfrow=c(1,1))
matplot(((a1$burn.in:a1$Tmax)-a1$burn.in)/52, t(a1$storeI[,a1$burn.in:a1$Tmax]),type="l", xlab="", ylab="cases", lty=1, col=cols) #,log="y"
matplot(((a2$burn.in:a2$Tmax)-a2$burn.in)/52, t(a2$storeI[,a1$burn.in:a1$Tmax]),type="l", add=TRUE,lty=3, col=cols) 
abline(v=1:20, lty=3, col="grey",lwd=0.5)

## we now need to convert cases into age.... might need age strcuture after all :) 




