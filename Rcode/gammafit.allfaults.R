# This file sources gammafit.R and fits a gamma renewal process 
# to the occurrence times from each fault segment.

## Set the Working Directory to the "Rcode" folder
#setwd("/PaleoEQ/Rcode")

# Load data
a <- dir("../DataFinal/chronologies_all_final")
nfault <- length(a)

K = 100

## Setting the values for n.inter, n.burnin, n.thin and gelman.cutoff for the jags run. 
## The values used in the manuscript are 
n.iter.mc = 5010000 
n.burnin.mc = 10000
n.thin.mc = 1000 
gelman.cutoff = 1.02
## which are the default values in the bptfit.R file. 
## If one wants to test if the code runs properly, then use the following three
## lines to reduce the computational time, which most likely will not result in
## convergence of the MCMC chains, but can test that the code is running.
#n.iter.mc = 6000
#n.burnin.mc = 1000
#n.thin.mc = 100
#gelman.cutoff = 2

## For convergence of MCMC chains, a minimum of the following values are suggested.
# n.iter.mc = 55000
# n.burnin.mc = 5000
# n.thin.mc = 10
# gelman.cutoff = 1.2


gelman <- matrix(NA,nrow=nfault,ncol=7)


############################################################################
#### Uncomment the following lines if one round of this file has been run
#### but some fault segments have the Gelman diagnostics values greater than 
#### 1.02, that is the chains are not converged for those segments.
#### Then repeat running this file until all chains are converged for 
#### all fault segments.
###########################################################################
## load Gelman Diagnostics and check which faults didn't converge

#b <- dir("../Results/GammaRes/GelmanDiag")
#bsub <- substr(b,16,17)
#bsub <- sort(as.numeric(bsub))

#gelmandiag <- matrix(NA,nrow=length(bsub),ncol=7)
#gelmanmax <- c()
#count <- 1
#for (j in bsub){
# temp <- read.csv(paste("../Results/GammaRes/GelmanDiag/GelmanDiagGamma",j,".csv",sep=""))
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

# i = as.integer(system('echo $KK', intern=T))


## If there is no access to supercomputers, then one can use a loop on a PC
## but this will take a long time to finish if n.iter = 5010000, n.burnin = 10000,
## n.thin=1000 are used in the gammafit.R file.

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

# Calculate inter-event times. 
# The last column is the censored time - last event occurrence time
# isCensor: each row contains repeated 0's for interevent times and 1 for
# censored times
inter.cens = t(diff(t(t.occ.cens))) # Transpose, take difference and transpose back
isCensor = t(replicate(m,c(rep(0,N-1),1)))


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


source('../Rcode/gammafit.R')


s1 = Sys.time()
gelman.1.2 <- FALSE
while(!gelman.1.2){
catched <- FALSE 
while(!catched) 
{
 tryCatch({
   mod <- gammafit(i,inter.NA, t.occ, N, K, isCensor, CensLim,
                   n.iter.mc,n.burnin.mc,n.thin.mc)
   catched = TRUE
   }, error=function(e){})
 print(catched)
}
gelman <- gelman.diag(mod[,c(1:6,6+K)])$psrf[,1]
if (all(gelman<=gelman.cutoff)) gelman.1.2 = TRUE
}


mod.mat <- as.matrix(mod)
mod.dat <- as.data.frame(mod.mat)
t.occ.Fore.vec <- as.vector(mod.mat[,6:(5+K)])


library("MCMCvis")
hpd <- MCMCsummary(mod,HPD=TRUE,hpd_prob=0.95,round=3)

source("../Rcode/hpd.interval.R")
hpd.t.occ.N <- hpd.interval(t.occ.Fore.vec,prob=0.95)

hpd.gamma <- list("hpd"=hpd,"hpd.t.occ.N"=hpd.t.occ.N)


t.fore.mean <- as.vector(mod.mat[,6+K])


write.csv(t(gelman),file=paste("../Results/GammaRes/GelmanDiag/GelmanDiagGamma",i,".csv",sep=""))
save(hpd.gamma,file=paste("../Results/GammaRes/HPDGamma",i,".image",sep=""))

}



q()


