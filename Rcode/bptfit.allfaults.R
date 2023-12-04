# This file sources bptfit.R and fits a bpt renewal process 
# to the occurrence times from each fault segment.

# Load data
a <- dir("../DataFinal/chronologies_all_final")
nfault <- length(a)

K = 100

## Setting the values for n.inter, n.burnin, n.thin and gelman.cutoff for the jags run. 
## The values used in the manuscript are 
# n.iter.mc = 5010000 
# n.burnin.mc = 10000
# n.thin.mc=1000 
# gelman.cutoff <- 1.02
## which are the default values in the bptfit.R file. 
## If one wants to test if the code runs properly, then use the following three
## lines to reduce the computational time, which most likely will not result in
## convergence of the MCMC chains, but can test that the code is running.
n.iter.mc = 6000
n.burnin.mc = 1000
n.thin.mc=1
gelman.cutoff <- 2

## For convergence of MCMC chains, a minimum of the following values are suggested.
# n.iter.mc = 55000
# n.burnin.mc = 5000
# n.thin.mc=10
# gelman.cutoff <- 1.2


gelman <- matrix(NA,nrow=nfault,ncol=7)


############################################################################
#### Uncomment the following lines if one round of this file has been run
#### but some fault segments have the Gelman diagnostics values greater than 
#### 1.02, that is the chains are not converged for those segments.
#### Then repeat running this file until all chains are converged for 
#### all fault segments.
###########################################################################
## load Gelman Diagnostics and check which faults didn't converge

#b <- dir("../Results/bptRes/GelmanDiag")
#bsub <- substr(b,14,15)
#bsub <- sort(as.numeric(bsub))

#gelmandiag <- matrix(NA,nrow=length(bsub),ncol=7)
#gelmanmax <- c()
#count <- 1
#for (j in bsub){
# temp <- read.csv(paste("../Results/bptRes/GelmanDiag/GelmanDiagBPT",j,".csv",sep=""))
# gelmandiag[count,] <- as.numeric(temp[-1]) 
# gelmanmax[count] <- max(as.numeric(temp[-1]))
# count <- count+1
#}

#if (length(which(gelmanmax>=1.02))>=1){
#  ind <- bsub[-which(gelmanmax>=1.02)]
#}else {ind <- bsub}

#index <- (1:93)[-ind]
#i = index[as.integer(system('echo $KK', intern=T))]
############################################################################
#### Uncomment the lines above if one round of this file has been run
#### but some fault segments have the Gelman diagnostics values greater than 
#### 1.02, that is the chains are not converged for those segments.
###########################################################################



## This line is used if we use supercomputer .sl file to run different fault 
## segments in parallel.

#i = as.integer(system('echo $KK', intern=T))




## If there is no access to supercomputers, then one can use a loop on a PC
## but this will take a long time to finish if n.iter.mc = 5010000, 
## n.burnin.mc = 10000, n.thin.mc=1000 are used.

for (i in 1:nfault){

data = read.csv(paste('../DataFinal/chronologies_all_final/',a[i],sep=''), header=FALSE)


# Chronologies have occurrence times in calendar years. 
# m: number of MC samples of the chronologies
# N: Number of earthquakes at the fault segment 
t.occ = data.matrix(data)
m = nrow(t.occ)
N <- ncol(t.occ)

# Add a column with 2022 to indicate
# the censored time at the year 2022, for forecasting purpose
t.occ.cens = cbind(t.occ,rep(2022,m))

## Calculate inter-event times. 
## The last column is the censored time - last event occurrence time
## isCensor: each row contains repeated 0's for interevent times and 1 for
## censored times
inter.cens = t(diff(t(t.occ.cens))) # Transpose, take difference and transpose back
isCensor = t(replicate(m,c(rep(0,N-1),1)))
max.inter <- max(inter.cens)

## inter.NA contains the interevent times and the last column is NA
## as the next earthquake hasn't occurred and we want to forecast its 
## occurrence time. JAGS will simulate the next occurrence time using MCMC
inter.NA = inter.cens
inter.NA[,N] = NA

  
# CensLim: KxN matrix. Censoring time for each element of inter.NA. 
### If the event was observed for inter.NA[i,j], then the "censoring
### time" must be greater than the observed inter-event time (inter.NA[i,j]). 
### We use inter-event times + 1 for non-censored data and use
### (year 2022 - last observed earthquake occurrence time for the ith MC sample)
### as censoring time for inter.NA[i,N]

isCensor10 = 1 - isCensor
CensLim = inter.cens + isCensor10


source('../Rcode/bptfit.R')


###########  Comment out the following two lines to not run MCMC again but just read the MCMC results
s1 = Sys.time()
gelman.1.2 <- FALSE
while(!gelman.1.2){
catched <- FALSE 
while(!catched) 
{
 tryCatch({
   mod <- bptfit(i,inter.NA, N, K, CensLim,n.iter.mc,n.burnin.mc,n.thin.mc)
   catched = TRUE
   }, error=function(e){})
 print(catched)
}
gelman <- gelman.diag(mod[,c(1:6,6+K)])$psrf[,1]
if (all(gelman<=gelman.cutoff)) gelman.1.2 = TRUE
}
###########


###########  Comment out the following four lines if running MCMC 
#library(lattice)
#library(R2jags)	
#load(paste("../Results/bptRes/bptFault-mcmc",i,".image",sep=""))
#mod <- mod.mcmc
###########


mod.mat <- as.matrix(mod)
mod.dat <- as.data.frame(mod.mat)


library(actuar)

## P(X <= x | X > 2022 - t_n)
Fx.cond <- function(q,mean0,shape0,t.cens){
  (pinvgauss(q, mean=mean0, shape = shape0)-
     pinvgauss(t.cens, mean=mean0, 
               shape = shape0))/
    pinvgauss(t.cens, mean=mean0, 
              shape = shape0,lower.tail=FALSE)
}

## Function to solve F(x|X>t.cens) = runif(1,0,1)
Finv <- function(mean0, shape0, t.cens, ran.num){
  ran.inv <- function(q, ran.num){
    Fx.cond(q,mean0,shape0,t.cens) - ran.num
  }
  uniroot(ran.inv, interval = c(t.cens,max.inter*10), ran.num = ran.num)
}


# Generate forecast occurrence times from the posterior samples
n.num <- length(mod.dat$mu)
inter.fore <- minval <- matrix(NA,n.num,K)
t.occ.N <- matrix(NA,n.num,K)


for (m in 1:K){
  Ym <- eval(parse(text=paste('mod.dat$"Y[',m,']"',sep="")))
  Zm <- eval(parse(text=paste('mod.dat$"Z[',m,']"',sep="")))
  mu.p <- mod.dat$mu*Zm
  shape.p <- mod.dat$mu*Zm/(mod.dat$alpha*Ym)^2
  unif.num <- runif(n.num,0,1)

  for (ii in 1:n.num){
    catched <- FALSE 
    while(!catched) 
    {
      ran.num <- runif(1, 0, 1)
      t.cens <- 2022-t.occ[m,N]

      tryCatch({
       ## Forecast inter-event times
        inter.fore[ii,m] <- Finv(mu.p[ii], shape.p[ii], t.cens, ran.num)$root 
        catched = TRUE
      }, error=function(e){})
      print(catched)
    }

    t.occ.N[ii,m] <- t.occ[m,N]+inter.fore[ii,m]
  }
  print(m)
}


t.occ.Fore.vec <- as.vector(t.occ.N)

s2 = Sys.time()
s2-s1


library("MCMCvis")
hpd <- MCMCsummary(mod,HPD=TRUE,hpd_prob=0.95,round=3)


source("../Rcode/hpd.interval.R")
hpd.t.occ.N <- hpd.interval(t.occ.Fore.vec,prob=0.95)

hpd.bpt <- list("hpd"=hpd,"hpd.t.occ.N"=hpd.t.occ.N)


}



bpt.res <- list("hpd.bpt"=hpd.bpt,"t.occ.N"=t.occ.N) 

write.csv(t(gelman),file=paste("../Results/bptRes/GelmanDiagBPT",i,".csv",sep=""))
save(bpt.res,file=paste("../Results/bptRes/HPD-toccBPT",i,".image",sep=""))


q()




