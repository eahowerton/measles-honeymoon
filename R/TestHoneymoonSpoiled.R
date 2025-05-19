library(epimdr2)

source("R/AgeStructModel-RSV.R")


 ## series from flat to peaky to diagonal to polymod
 background <- 0.001
  waifw1 <- matrix(1,length(c(1:120,seq(180,600,by=60))/12),length(c(1:120,seq(180,600,by=60))/12))
  waifw2 <- get.smooth.WAIFW(age.class.boundries=c(1:120,seq(180,600,by=60))/12, mu=5,sig=0.2, gam=0.05, delta= background)
  waifw3 <- get.smooth.WAIFW(age.class.boundries=c(1:120,seq(180,600,by=60))/12, mu=10,sig=0.1, gam=0.01, delta= background)
  waifw4 <- get.smooth.WAIFW(age.class.boundries=c(1:120,seq(180,600,by=60))/12, mu=7,sig=0.6, gam=0.05, delta= background) 
  waifw5 <- get.polymod.WAIFW(age.class.boundries=c(1:120,seq(180,600,by=60))/12, do.touch=TRUE)


  ## time controllers
  n.per.year <- 26
  time.step <- 1/2
 
  ## basic demography
  age.classes <- c(1:120,seq(180,600,by=60))
  mort <- c(rep(0,127),1)
  fert <-  c(rep(0,127),1) ##note that if lambda is not 1; vast slippage in average age
  pop.n <- 500000
  waning.maternal <- 0.99 #very rapid TODO CHECK
    
  ## tweak this outside of this function to adjust the rate of increase to get something close to equilbirium (corresponds to lambda=1)
  ## here, using this to get reasonable starting values
  xx <- findStableStruct(age.classes=age.classes, 
                         mort=mort,fert =fert,time.step = time.step)
  #xx$lambda



  ## pick out R0 and the desired WAIFW (available R0 for scaling is 10 to 24)
  R0 <- 14
  chs.waifw <- 2
  waifw.adj <- waifw2

  ## make the waifw what it 'should' be for this R0 (this step is used in aFullComparisonSimple so scalar will depend on it
  waifw.adj <- scale.Waifw(R0=R0, DFE=xx$stable.age*pop.n, waifw=waifw.adj)

  ##bring in the multpliers: 
  find.multipliers<-read.csv("~/Library/CloudStorage/OneDrive-PrincetonUniversity/Documents/temp/adj.valsFALSE.csv",row.names=NULL)
  find.multipliers <- find.multipliers[find.multipliers[,1]==1,] ## pick out just the incidence adjusted. 
  find.multipliers <- find.multipliers[find.multipliers[,3]==R0,] ## pick out focal R0
  find.multipliers <- find.multipliers[find.multipliers[,2]==chs.waifw,] ## pick out focal waifw
  print(as.numeric(find.multipliers$V3))


  ## adjust the waifw
  waifw.adj <- waifw.adj*as.numeric(find.multipliers$V3)

  ### If all the vaccination happens in one month (two bisteps), and the target is 0.8, then the proportion vaccinated should be: log(1-0.8)/(-2)
  ### Note that indivdiuals are also aging out, which complicates this a bit. Effectively 3/4 stay in the age class, so if cover = 0.9 actually getting just get 0.67
  ### For this reason, reasonable to vaccinate in two montly age classes....


## DO some trials ####################

vacc1 <- 0.5
vacc2 <- 0.1

 ## Keep vaccination the same 
   tmp <- iterateSIRvacc(age.classes=age.classes, mort=mort,	
                    waifw=waifw.adj,
		    vacc.cover.start= vacc1, vacc.cover.release= vacc1,
		    alpha=0.6, time.step=time.step,
                    fert= fert,
                    burnin=10* n.per.year,
		    Tmax=(150* n.per.year), 
                    do.plot=F,start.pop=c(),start.pop.size=pop.n,
		    start.full.pop.struct=c(),gamma=0.97, do.unity=T)



  ## release in the honeymoon
   tmp1 <- iterateSIRvacc(age.classes=age.classes, mort=mort,	
                    waifw=waifw.adj,
		    vacc.cover.start= vacc1, vacc.cover.release= vacc2,
		    alpha=0.6, time.step=time.step,
                    fert= fert,
                    burnin=10* n.per.year,
		    Tmax=(150* n.per.year), 
                    do.plot=F,start.pop=c(),start.pop.size=pop.n,
		    start.full.pop.struct=c(),gamma=0.97, do.unity=T)


  age.vals <- extractAge(tmp$nt.store[tmp$i.inds,],tmp$age.classes)
  age.vals1 <- extractAge(tmp1$nt.store[tmp$i.inds,],tmp$age.classes)
  

###
chs.times <- 2104:2450
all.seas <- rep(tmp$seas,ncol(tmp$nt.store)/26)

par(mfrow=c(3,1), bty="l",pty="m", mar=c(4,5,1,2))
 plot(colSums(tmp$nt.store[tmp$i.inds, chs.times]), type="l",xlab="", ylab="Infected",lwd=2, cex.lab=1.5)#,log="y")
 abline(v=26*c(1:1000),lty=3)
 abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times]!=tmp1$vacc.vals[-1][chs.times])[1], col="red",lwd=3)

 plot(colSums(tmp$nt.store[tmp$s.inds, chs.times]), type="l",xlab="", ylab="Susceptible",lwd=2, cex.lab=1.5)#,log="y")
 abline(v=26*c(1:1000),lty=3)
 abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times]!=tmp1$vacc.vals[-1][chs.times])[1], col="red",lwd=3)

  ## unity beta but now try and pull away the seasonal swing
  plot(tmp$store.beta.unity[chs.times]/all.seas[chs.times], type="l",xlab="", ylab="Unity beta",lwd=2, cex.lab=1.5)#,log="y")
  points(tmp$store.beta.unity[chs.times], type="l",col="black",lty=4)
  points(all.seas[chs.times]*max(tmp$store.beta.unity[chs.times])/3,type="l",col="grey")
  abline(v=26*c(1:1000),lty=3)
  abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times]!=tmp1$vacc.vals[-1][chs.times])[1], col="red",lwd=3)



par(mfrow=c(3,1), bty="l",pty="m", mar=c(4,5,1,2))
 plot(colSums(tmp$nt.store[tmp$i.inds, chs.times]), type="l",xlab="", ylab="Infected",lwd=2, cex.lab=1.5)#,log="y")
 points(colSums(tmp1$nt.store[tmp$i.inds, chs.times]), type="l",col="blue")
 abline(v=26*c(1:1000),lty=3)
 abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times]!=tmp1$vacc.vals[-1][chs.times]), col="red",lwd=3)

 plot(colSums(tmp$nt.store[tmp$s.inds, chs.times]), type="l",xlab="", ylab="Susceptible",lwd=2, cex.lab=1.5)#,log="y")
  points(colSums(tmp1$nt.store[tmp$s.inds, chs.times]), type="l",col="blue")
 abline(v=26*c(1:1000),lty=3)
 abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times]!=tmp1$vacc.vals[-1][chs.times]), col="red",lwd=3)

  ## unity beta but now try and pull away the seasonal swing
  plot(tmp$store.beta.unity[chs.times]/all.seas[chs.times], type="l",xlab="", ylab="Unity beta",lwd=2, cex.lab=1.5)#,log="y")
  points(tmp1$store.beta.unity[chs.times]/all.seas[chs.times], type="l",col="blue")
  points(tmp$store.beta.unity[chs.times], type="l",col="black",lty=4)
  points(tmp1$store.beta.unity[chs.times], type="l",col="blue",lty=4)
  points(all.seas[chs.times]*max(tmp$store.beta.unity[chs.times])/3,type="l",col="grey")
  abline(v=26*c(1:1000),lty=3)
  abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times]!=tmp1$vacc.vals[-1][chs.times]), col="red",lwd=3)






chs.times1 <- chs.times
chs.times1 <- 2104: 2800

par(mfrow=c(3,1), bty="l",pty="m", mar=c(4,5,1,2))
 plot(colSums(tmp$nt.store[tmp$i.inds, chs.times1]), type="l",xlab="", ylab="Infected",lwd=2, cex.lab=1.5)#,log="y")
 points(colSums(tmp1$nt.store[tmp$i.inds, chs.times1]), type="l",col="blue")
 abline(v=26*c(1:1000),lty=3)
 abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times1]!=tmp1$vacc.vals[-1][chs.times1]), col="red",lwd=3)

 plot(tmp$store.beta.unity[chs.times1]/all.seas[chs.times1], type="l",xlab="", ylab="Unity beta",lwd=2, cex.lab=1.5)#,log="y")
  points(tmp1$store.beta.unity[chs.times1]/all.seas[chs.times1], type="l",col="blue")
  #points(tmp$store.beta.unity[chs.times1], type="l",col="black",lty=4)
  #points(tmp1$store.beta.unity[chs.times1], type="l",col="blue",lty=4)
  points(all.seas[chs.times1]*max(tmp$store.beta.unity[chs.times1])/3,type="l",col="grey")
  abline(v=26*c(1:1000),lty=3)
  abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times1]!=tmp1$vacc.vals[-1][chs.times1]), col="red",lwd=3)

plot(age.vals$avg.age[chs.times1]/12, type="l",xlab="", ylab="Age",lwd=2, cex.lab=1.5)#,log="y")
 points(age.vals1$avg.age[chs.times1]/12, type="l",col=4)
 abline(v=26*c(1:1000),lty=3)
 abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times1]!=tmp1$vacc.vals[-1][chs.times1]), col="red",lwd=3)














  ## pick out R0 and the desired WAIFW (available R0 for scaling is 10 to 24)
  R0 <- 14
  chs.waifw <- 2
  waifw.adj <- waifw2

  ## make the waifw what it 'should' be for this R0 (this step is used in aFullComparisonSimple so scalar will depend on it
  waifw.adj <- scale.Waifw(R0=R0, DFE=xx$stable.age*pop.n, waifw=waifw.adj)

  ##bring in the multpliers: 
  find.multipliers<-read.csv("~/Library/CloudStorage/OneDrive-PrincetonUniversity/Documents/temp/adj.valsFALSE.csv",row.names=NULL)
  find.multipliers <- find.multipliers[find.multipliers[,1]==1,] ## pick out just the incidence adjusted. 
  find.multipliers <- find.multipliers[find.multipliers[,3]==R0,] ## pick out focal R0
  find.multipliers <- find.multipliers[find.multipliers[,2]==chs.waifw,] ## pick out focal waifw
  print(as.numeric(find.multipliers$V3))


  ## adjust the waifw
  waifw.adj <- waifw.adj*as.numeric(find.multipliers$V3)

  ### If all the vaccination happens in one month (two bisteps), and the target is 0.8, then the proportion vaccinated should be: log(1-0.8)/(-2)
  ### Note that indivdiuals are also aging out, which complicates this a bit. Effectively 3/4 stay in the age class, so if cover = 0.9 actually getting just get 0.67
  ### For this reason, reasonable to vaccinate in two montly age classes....


## DO some trials ####################

vacc1 <- 0.5
vacc2 <- 0.2

 ## Keep vaccination the same 
   tmp <- iterateSIRvacc(age.classes=age.classes, mort=mort,	
                    waifw=waifw.adj,
		    vacc.cover.start= vacc1, vacc.cover.release= vacc1,
		    alpha=0.6, time.step=time.step,
                    fert= fert,
                    burnin=10* n.per.year,
		    Tmax=(150* n.per.year), 
                    do.plot=F,start.pop=c(),start.pop.size=pop.n,
		    start.full.pop.struct=c(),gamma=0.97, do.unity=T)



  ## release in the honeymoon
   tmp1 <- iterateSIRvacc(age.classes=age.classes, mort=mort,	
                    waifw=waifw.adj,
		    vacc.cover.start= vacc1, vacc.cover.release= vacc2,
		    alpha=0.6, time.step=time.step,
                    fert= fert,
                    burnin=10* n.per.year,
		    Tmax=(150* n.per.year), 
                    do.plot=F,start.pop=c(),start.pop.size=pop.n,
		    start.full.pop.struct=c(),gamma=0.97, do.unity=T)


  age.vals <- extractAge(tmp$nt.store[tmp$i.inds,],tmp$age.classes)
  age.vals1 <- extractAge(tmp1$nt.store[tmp$i.inds,],tmp$age.classes)
  

###
chs.times <- 2104:2450
all.seas <- rep(tmp$seas,ncol(tmp$nt.store)/26)

par(mfrow=c(3,1), bty="l",pty="m", mar=c(4,5,1,2))
 plot(colSums(tmp$nt.store[tmp$i.inds, chs.times]), type="l",xlab="", ylab="Infected",lwd=2, cex.lab=1.5)#,log="y")
 abline(v=26*c(1:1000),lty=3)
 abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times]!=tmp1$vacc.vals[-1][chs.times])[1], col="red",lwd=3)

 plot(colSums(tmp$nt.store[tmp$s.inds, chs.times]), type="l",xlab="", ylab="Susceptible",lwd=2, cex.lab=1.5)#,log="y")
 abline(v=26*c(1:1000),lty=3)
 abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times]!=tmp1$vacc.vals[-1][chs.times])[1], col="red",lwd=3)

  ## unity beta but now try and pull away the seasonal swing
  plot(tmp$store.beta.unity[chs.times]/all.seas[chs.times], type="l",xlab="", ylab="Unity beta",lwd=2, cex.lab=1.5)#,log="y")
  points(tmp$store.beta.unity[chs.times], type="l",col="black",lty=4)
  points(all.seas[chs.times]*max(tmp$store.beta.unity[chs.times])/3,type="l",col="grey")
  abline(v=26*c(1:1000),lty=3)
  abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times]!=tmp1$vacc.vals[-1][chs.times])[1], col="red",lwd=3)



par(mfrow=c(3,1), bty="l",pty="m", mar=c(4,5,1,2))
 plot(colSums(tmp$nt.store[tmp$i.inds, chs.times]), type="l",xlab="", ylab="Infected",lwd=2, cex.lab=1.5)#,log="y")
 points(colSums(tmp1$nt.store[tmp$i.inds, chs.times]), type="l",col="blue")
 abline(v=26*c(1:1000),lty=3)
 abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times]!=tmp1$vacc.vals[-1][chs.times]), col="red",lwd=3)

 plot(colSums(tmp$nt.store[tmp$s.inds, chs.times]), type="l",xlab="", ylab="Susceptible",lwd=2, cex.lab=1.5)#,log="y")
  points(colSums(tmp1$nt.store[tmp$s.inds, chs.times]), type="l",col="blue")
 abline(v=26*c(1:1000),lty=3)
 abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times]!=tmp1$vacc.vals[-1][chs.times]), col="red",lwd=3)

  ## unity beta but now try and pull away the seasonal swing
  plot(tmp$store.beta.unity[chs.times]/all.seas[chs.times], type="l",xlab="", ylab="Unity beta",lwd=2, cex.lab=1.5)#,log="y")
  points(tmp1$store.beta.unity[chs.times]/all.seas[chs.times], type="l",col="blue")
  points(tmp$store.beta.unity[chs.times], type="l",col="black",lty=4)
  points(tmp1$store.beta.unity[chs.times], type="l",col="blue",lty=4)
  points(all.seas[chs.times]*max(tmp$store.beta.unity[chs.times])/3,type="l",col="grey")
  abline(v=26*c(1:1000),lty=3)
  abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times]!=tmp1$vacc.vals[-1][chs.times]), col="red",lwd=3)






chs.times1 <- chs.times
chs.times1 <- 2104: 2800

par(mfrow=c(3,1), bty="l",pty="m", mar=c(4,5,1,2))
 plot(colSums(tmp$nt.store[tmp$i.inds, chs.times1]), type="l",xlab="", ylab="Infected",lwd=2, cex.lab=1.5)#,log="y")
 points(colSums(tmp1$nt.store[tmp$i.inds, chs.times1]), type="l",col="blue")
 abline(v=26*c(1:1000),lty=3)
 abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times1]!=tmp1$vacc.vals[-1][chs.times1]), col="red",lwd=3)

 plot(tmp$store.beta.unity[chs.times1]/all.seas[chs.times1], type="l",xlab="", ylab="Unity beta",lwd=2, cex.lab=1.5)#,log="y")
  points(tmp1$store.beta.unity[chs.times1]/all.seas[chs.times1], type="l",col="blue")
  #points(tmp$store.beta.unity[chs.times1], type="l",col="black",lty=4)
  #points(tmp1$store.beta.unity[chs.times1], type="l",col="blue",lty=4)
  points(all.seas[chs.times1]*max(tmp$store.beta.unity[chs.times1])/3,type="l",col="grey")
  abline(v=26*c(1:1000),lty=3)
  abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times1]!=tmp1$vacc.vals[-1][chs.times1]), col="red",lwd=3)

plot(age.vals$avg.age[chs.times1]/12, type="l",xlab="", ylab="Age",lwd=2, cex.lab=1.5)#,log="y")
 points(age.vals1$avg.age[chs.times1]/12, type="l",col=4)
 abline(v=26*c(1:1000),lty=3)
 abline(v=which(tmp1$vacc.vals[-length(tmp1$vacc.vals)][chs.times1]!=tmp1$vacc.vals[-1][chs.times1]), col="red",lwd=3)


