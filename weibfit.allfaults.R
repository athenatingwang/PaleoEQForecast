# This file sources weibfit.R and fits
# a Weibull renewal process to the occurrence times from each fault.

# Load data
a <- dir("chronologies_all_final")
nfault <- length(a)

K = 100
#gelman <- matrix(NA,nrow=nfault,ncol=7)

## load Gelman Diagnostics and check which faults didn't converge
b <- dir("/nesi/nobackup/nesi00326/WeibRes/GelmanDiag")
bsub <- substr(b,15,16)
bsub <- sort(as.numeric(bsub))

gelmandiag <- matrix(NA,nrow=length(bsub),ncol=7)
gelmanmax <- c()
count <- 1
for (j in bsub){
 temp <- read.csv(paste("/nesi/nobackup/nesi00326/WeibRes/GelmanDiag/GelmanDiagWeib",j,".csv",sep=""))
 gelmandiag[count,] <- as.numeric(temp[-1]) 
 gelmanmax[count] <- max(as.numeric(temp[-1]))
 count <- count+1
}

if (length(which(gelmanmax>=1.02))>=1){
  ind <- bsub[-which(gelmanmax>=1.02)]
}else {ind <- bsub}

index <- (1:93)[-ind]



i = index[as.integer(system('echo $KK', intern=T))]



#i = as.integer(system('echo $KK', intern=T))


#for (i in 1:nfault){

data = read.csv(paste('chronologies_all_final/',a[i],sep=''), header=FALSE)


#pdf(paste('Weib',i,'.pdf',sep=''))

# Chronologies have occurrence times in calendar years. 
# m: number of MC samples of the chronologies
# N: Number of earthquakes at the fault
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

# CensLim: Element i of each row is the censoring time for case i. 
# If the event was observed for case i then the "censoring
# time" must be greater than the observed event time. 
# We use inter-event times + 1 for non-censored data and 
# the censored time from the last censored time

#isCensor10 = ifelse(isCensor==1,0,1)
isCensor10 = 1 - isCensor
CensLim = inter.cens + isCensor10

# inter.NA contains the interevent times and the last column is NA
# as the next earthquake hasn't occurred and we want to forecast its 
# occurrence time. JAGS will simulate the next occurrence time using MCMC
inter.NA = inter.cens
inter.NA[,N] = NA
source('weibfit.R')



###########  Comment out the following two lines to not run MCMC again but just read the MCMC results
gelman.1.2 <- FALSE
while(!gelman.1.2){
catched <- FALSE 
while(!catched) 
{
 tryCatch({
   mod <- weibfit(i,inter.NA, t.occ, N, K, isCensor, CensLim, plot.Y=T)
   catched = TRUE
   }, error=function(e){})
 print(catched)
}
gelman <- gelman.diag(mod[,c(1:6,6+K)])$psrf[,1]
if (all(gelman<=1.02)) gelman.1.2 = TRUE
}
###########




###########  Comment out the following four lines if running MCMC 
#library(lattice)
#library(R2jags)	
#load(paste("/nesi/nobackup/nesi00326/WeibRes/WeibFault-mcmc",i,".image",sep=""))
#mod <- mod.mcmc
###########



mod.mat <- as.matrix(mod)
mod.dat <- as.data.frame(mod.mat)
t.occ.Fore.vec <- as.vector(mod.mat[,6:(5+K)])


library("MCMCvis")
hpd <- MCMCsummary(mod,HPD=TRUE,hpd_prob=0.95,round=3)

source("hpd.interval.R")
hpd.t.occ.N <- hpd.interval(t.occ.Fore.vec,prob=0.95)

hpd.weib <- list("hpd"=hpd,"hpd.t.occ.N"=hpd.t.occ.N)


# Some more plots
#xyplot(mod[,c(1:6,6+K)], layout=c(4,2), aspect="fill")

# Density plot
#densityplot(mod[,c(1:6,6+K)], layout=c(4,2), aspect="fill")


#Auto-correlation plot
#autocorr.plot(mod[,c(1:6,6+K)])


x.lim <- c(min(t.occ.Fore.vec),max(t.occ.Fore.vec))
#par(mfrow=c(2,1))
#dens <- density(t.occ.Fore.vec)
#plot(dens,main="Density of forecast interevent times",
#     xlim=x.lim)
#lines(t.occ.Fore.vec, rep(max(dens$y)/100, length(t.occ.Fore.vec)), 
#      type = "h")

t.fore.mean <- as.vector(mod.mat[,6+K])
#dens <- density(t.fore.mean)
#plot(dens,main="Density of mean forecast interevent times",
#     xlim=x.lim)
#lines(t.fore.mean, rep(max(dens$y)/100, length(t.fore.mean)), 
#      type = "h")

#dev.off()

#}



write.csv(t(gelman),file=paste("/nesi/nobackup/nesi00326/WeibRes/GelmanDiagWeib",i,".csv",sep=""))
save(hpd.weib,file=paste("/nesi/nobackup/nesi00326/WeibRes/HPDWeib",i,".image",sep=""))


q()



