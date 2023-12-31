# This file sources poisfitNoCensor.R and fits a Poisson process 
# to the occurrence times from each fault segment.

## Set the Working Directory to the "/Demo/RetroForecast" folder
#setwd("/PaleoEQ/Demo/RetroForecast")

# Load data
a <- dir("../../DataFinal/chronologies_all_final")
nfault <- length(a)

K = 50

## Setting the values for n.inter, n.burnin, n.thin and gelman.cutoff for the jags run. 
## The values used in the manuscript are 
#n.iter.mc = 5010000 
#n.burnin.mc = 10000
#n.thin.mc=1000 
#gelman.cutoff <- 1.02
## which are the default values in the bptfit.R file. 
## If one wants to test if the code runs properly, then use the following three
## lines to reduce the computational time, which most likely will not result in
## convergence of the MCMC chains, but can test that the code is running.
n.iter.mc = 6000
n.burnin.mc = 1000
n.thin.mc = 100
gelman.cutoff = 5

## For convergence of MCMC chains, a minimum of the following values are suggested.
# n.iter.mc = 55000
# n.burnin.mc = 5000
# n.thin.mc=10
# gelman.cutoff <- 1.2


gelman <- matrix(NA,nrow=nfault,ncol=5)


############################################################################
#### Uncomment the following lines if one round of this file has been run
#### but some fault segments have the Gelman diagnostics values greater than 
#### 1.02, that is the chains are not converged for those segments.
#### Then repeat running this file until all chains are converged for 
#### all fault segments.
###########################################################################
## load Gelman Diagnostics and check which faults didn't converge

#b <- dir("../../Results/PoisSimRes2/GelmanDiag")
#bsub <- substr(b,15,16)
#bsub <- sort(as.numeric(bsub))

#gelmandiag <- matrix(NA,nrow=length(bsub),ncol=5)
#gelmanmax <- c()
#count <- 1
#for (j in bsub){
# temp <- read.csv(paste("../../Results/PoisSimRes2/GelmanDiag/GelmanDiagPois",j,".csv",sep=""))
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
## but this will take a long time to finish if n.iter = 5010000, n.burnin = 10000,
## n.thin=1000 are used in the poisfitNoCensor.R file.

s1 = Sys.time()


for (i in 1:nfault){

data = read.csv(paste('../../DataFinal/chronologies_all_final/',a[i],sep=''), header=FALSE)


# Chronologies have occurrence times in calendar years. 
# m: number of MC samples of the chronologies
# N: Number of earthquakes at the fault
t.occ = data.matrix(data)
m = nrow(t.occ)
N1 <- ncol(t.occ)

# remove the last occurrence time and do NOT use censored time, 
# use the last occurrence time directly to forecast future values
# That is, we don't have a cutoff time
N <- N1-1
t.occ.cens = t.occ[,1:N]


# Calculate inter-event times. 
inter.cens = t(diff(t(t.occ.cens))) # Transpose, take difference and transpose back
  
# inter.NA contains the interevent times and the last column is NA
# as the next earthquake hasn't occurred and we want to forecast its 
# occurrence time. JAGS will simulate the next occurrence time using MCMC
inter.NA = matrix(NA,nrow=nrow(t.occ.cens),ncol=ncol(t.occ.cens))
inter.NA[,1:(N-1)] = inter.cens


source('../RetroForecast/poisfitNoCensor.R')


gelman.1.2 <- FALSE
while(!gelman.1.2){
catched <- FALSE 
while(!catched) 
{
 tryCatch({
   mod <- poisfit(i,inter.NA, t.occ, N, K,n.iter.mc,n.burnin.mc,n.thin.mc)
   catched = TRUE
   }, error=function(e){})
 print(catched)
}
gelman <- gelman.diag(mod[,c(1:4,4+K)])$psrf[,1]
if (all(gelman<=gelman.cutoff)) gelman.1.2 = TRUE
}


mod.mat <- as.matrix(mod)
mod.dat <- as.data.frame(mod.mat)
t.occ.Fore.vec <- as.vector(mod.mat[,4:(3+K)])


library("MCMCvis")
hpd <- MCMCsummary(mod,HPD=TRUE,hpd_prob=0.95,round=3)

source("../../Rcode/hpd.interval.R")
hpd.t.occ.N <- hpd.interval(t.occ.Fore.vec,prob=0.95)

hpd.pois <- list("hpd"=hpd,"hpd.t.occ.N"=hpd.t.occ.N)



write.csv(t(gelman),file=paste("../../ResultsDemo/PoisSimRes2/GelmanDiag/GelmanDiagPois",i,".csv",sep=""))
save(hpd.pois,file=paste("../../ResultsDemo/PoisSimRes2/HPDPois",i,".image",sep=""))


}

s2 = Sys.time()
s2-s1

#Time difference of 2.729989 mins

q()



