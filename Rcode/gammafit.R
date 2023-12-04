# Fitting a gamma renewal process to each MC sample of the earthquake chronologies 

# This file contains the jags model and general jags run for any data.
# Refer to gammafit.allfaults.R which sources this gammafit.R and fits
# the gamma renewal process to the occurrence times from each fault.


###############################################################
############               Jags model              ############
###############################################################
gamma.jags <- "
model{
  for (j in 1:K){ # jth MC sample of the data, K samples in total
   for (i in 1:N){ # ith observation of the jth MC sample	
     # Likelihood 
     isCensor[j,i] ~ dinterval( inter.NA[j,i] , CensLim[j,i] )
     inter.NA[j,i] ~ dgamma(Z[j]*alpha,Y[j]*beta)  
   }

   Z[j] ~ dgamma(sh, rt)
   Y[j] ~ dgamma(shy, rty)

   ## forecast the next occurrence time t.occ.N for the jth MC sample
   t.occ.N[j] <- t.occ[j,N] + inter.NA[j,N]  
  }

  beta ~ dnorm(0,10^-4)T(0,)
  alpha ~ dnorm(0,10^-4)T(0,)
  sh <- 1/pow(sigmaZ,2)
  rt <- 1/pow(sigmaZ,2)
  sigmaZ ~ dt(0,0.04,3)T(0,)
  shy <- 1/pow(sigmaY,2)
  rty <- 1/pow(sigmaY,2)
  sigmaY ~ dt(0,0.04,3)T(0,)

  # t.occ.N is the forecast time of the next earthquake occurrence
  tfore.mean <- mean(t.occ.N)
}
"
###############################################################
############               Jags model              ############
###############################################################




###################################################################################
#### Refer to gammafit.allfaults.R for how to get each input argument below from the 
#### earthquake chronologies.
###################################################################################
# faulti: integer indicating the ith fault in the folder chronologies_all_final
# N: Number of earthquakes at the fault
# K: Number of MC samples used from the earthquake chronologies (we used 100)

# inter.NA: KxN matrix. It contains the interevent times from the K MC samples, 
### and the last column is NA as the next earthquake hasn't occurred and we
### want to forecast its occurrence time. JAGS will simulate the next occurrence
### time using MCMC

# CensLim: KxN matrix. Censoring time for each element of inter.NA. 
### If the event was observed for inter.NA[i,j], then the "censoring
### time" must be greater than the observed inter-event time (inter.NA[i,j]). 
### We use inter-event times + 1 for non-censored data and use
### (year 2022 - last observed earthquake occurrence time for the ith MC sample)
### as censoring time for inter.NA[i,N]


gammafit <- function(faulti,inter.NA, t.occ, N, K, isCensor, CensLim,
                     n.iter.mc = 5010000,n.burnin.mc = 10000,n.thin.mc=1000){
	library(lattice)
	library(R2jags)	
	jagsdata <- list("inter.NA", "t.occ", "N", "K","isCensor","CensLim")
	
	# Define the parameters that we are interested in. 
	# Their posterior distributions will be the output.
	# t.occ.N is the forecast time of the next earthquake occurrence
  params <- c("alpha","beta","sigmaY","sigmaZ","t.occ.N","tfore.mean","Y","Z")

	# Set the initial values
	rate0 <- mean(inter.NA[,-N])/(mean(inter.NA[,-N]^2)-
	                       mean(inter.NA[,-N])^2)
	shape0 <- rate0*mean(inter.NA[,-N])
	
	inits <- function(){
			list("alpha"=shape0+runif(1,-shape0/4,shape0/4),"beta"=rate0+runif(1,-rate0/3,rate0/3))
			}

	mod <- jags(data = jagsdata, inits = inits,
		            parameters.to.save = params, n.chains = 3,
	              n.iter = n.iter.mc, n.burnin = n.burnin.mc,n.thin=n.thin.mc,
			    model.file = textConnection(gamma.jags))
	
	# Convert to an MCMC object
	mod.mcmc <- as.mcmc(mod)

### Can use the code below to save the result in a .image file		
#	save(mod.mcmc,file=paste("../Results/GammaRes/GammaFault-mcmc",faulti,".image",sep=""))

	return(mod.mcmc)
	}
	
	







