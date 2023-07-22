# This file sources bptfit.R and fits
# a bpt renewal process to the occurrence times from each fault.

# Load data
a <- dir("chronologies_all_final")
nfault <- length(a)

K = 100
#gelman <- matrix(NA,nrow=nfault,ncol=7)


## load Gelman Diagnostics and check which faults didn't converge

b <- dir("/nesi/nobackup/nesi00326/bptRes/GelmanDiag")
bsub <- substr(b,14,15)
bsub <- sort(as.numeric(bsub))

gelmandiag <- matrix(NA,nrow=length(bsub),ncol=7)
gelmanmax <- c()
count <- 1
for (j in bsub){
 temp <- read.csv(paste("/nesi/nobackup/nesi00326/bptRes/GelmanDiag/GelmanDiagBPT",j,".csv",sep=""))
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


#pdf(paste('bpt',i,'.pdf',sep=''))

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
max.inter <- max(inter.cens)
  
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
source('bptfit.R')


###########  Comment out the following two lines to not run MCMC again but just read the MCMC results
s1 = Sys.time()
gelman.1.2 <- FALSE
while(!gelman.1.2){
catched <- FALSE 
while(!catched) 
{
 tryCatch({
   mod <- bptfit(i,inter.NA, N, K, CensLim, plot.Y=T)
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
#load(paste("/nesi/nobackup/nesi00326/bptRes/bptFault-mcmc",i,".image",sep=""))
#mod <- mod.mcmc
###########

####mod <- bptfit(i,inter.NA, t.occ, N, K, isCensor, CensLim, plot.Y=T)



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
#    (Fx.cond(q,mean0,shape0,t.cens) - ran.num)^2
  }
  uniroot(ran.inv, interval = c(t.cens,max.inter*10), ran.num = ran.num)
#  nlm(ran.inv, runif(1,t.cens,max.inter*10), ran.num = ran.num)
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
  inter.fore[ii,m] <- Finv(mu.p[ii], shape.p[ii], t.cens, ran.num)$root ## Forecast inter-event times
catched = TRUE
   }, error=function(e){})
 print(catched)
}


#  inter.fore[ii,m] <- Finv(mu.p[ii], shape.p[ii], t.cens, ran.num)$estimate
#  minval[ii,m] <- Finv(mu.p[ii], shape.p[ii], t.cens, ran.num)$minimum
  t.occ.N[ii,m] <- t.occ[m,N]+inter.fore[ii,m]
}

print(m)
}



t.occ.Fore.vec <- as.vector(t.occ.N)

s2 = Sys.time()
s2-s1



library("MCMCvis")
hpd <- MCMCsummary(mod,HPD=TRUE,hpd_prob=0.95,round=3)


source("hpd.interval.R")
hpd.t.occ.N <- hpd.interval(t.occ.Fore.vec,prob=0.95)

hpd.bpt <- list("hpd"=hpd,"hpd.t.occ.N"=hpd.t.occ.N)


# Some more plots
#xyplot(mod[,c(1:6,6+K)], layout=c(4,2), aspect="fill")

# Density plot
#densityplot(mod[,c(1:6,6+K)], layout=c(4,2), aspect="fill")


#Auto-correlation plot
#autocorr.plot(mod[,c(1:6,6+K)])


x.lim <- c(min(t.occ.Fore.vec),max(t.occ.Fore.vec))
#par(mfrow=c(2,1))
#dens <- density(t.occ.Fore.vec)
#plot(dens,main="Density of forecast occurrence times",
#     xlim=x.lim)
#lines(t.occ.Fore.vec, rep(max(dens$y)/100, length(t.occ.Fore.vec)), 
#      type = "h")

#t.fore.mean <- as.vector(mod.mat[,6+K])
#dens <- density(t.fore.mean)
#plot(dens,main="Density of mean forecast occurrence times",
#     xlim=x.lim)
#lines(t.fore.mean, rep(max(dens$y)/100, length(t.fore.mean)), 
#      type = "h")

#dev.off()

#}



bpt.res <- list("hpd.bpt"=hpd.bpt,"t.occ.N"=t.occ.N) #,"minval"=minval)

write.csv(t(gelman),file=paste("/nesi/nobackup/nesi00326/bptRes/GelmanDiagBPT",i,".csv",sep=""))
save(bpt.res,file=paste("/nesi/nobackup/nesi00326/bptRes/HPD-toccBPT",i,".image",sep=""))


q()




