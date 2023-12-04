# Fitting a Poisson process without censoring to each MC sample of the earthquake chronologies 

# This file contains the jags model and general jags run for any data.
# Refer to poisfit.allfaults.R which sources this poisfitNoCensor.R and fits
# the Poisson process without censoring to the occurrence times from each fault segment.


###############################################################
############               Jags model              ############
###############################################################
poisNoCensor.jags <- "
model{
  for (j in 1:K){ # jth MC sample of the data, K samples in total
   for (i in 1:N){ # ith observation of the jth MC sample	
     # Likelihood 
     inter.NA[j,i] ~ dexp(lambda*Z[j]) 
   }
   Z[j] ~ dgamma(sh, rt)	

   ## forecast the next occurrence time t.occ.N for the jth MC sample
   t.occ.N[j] <- t.occ[j,N] + inter.NA[j,N]  
}
  sh <- 1/pow(sigmaZ,2)
  rt <- 1/pow(sigmaZ,2)
  sigmaZ ~ dt(0,0.04,3)T(0,)
  lambda ~ dnorm(0,10^-4)T(0,)

  # t.occ.N is the forecast time of the next earthquake occurrence
  tfore.mean <- mean(t.occ.N)
}
"
###############################################################
############               Jags model              ############
###############################################################




###################################################################################
#### Refer to poisfit.allfaults.R for how to get each input argument below from the 
#### earthquake chronologies.
###################################################################################
# faulti: integer indicating the ith fault in the folder chronologies_all_final
# N: Number of earthquakes at the fault
# K: Number of MC samples used from the earthquake chronologies (we used 100)

# inter.NA: KxN matrix. It contains the interevent times from the K MC samples, 
### and the last column is NA as the next earthquake hasn't occurred and we
### want to forecast its occurrence time. JAGS will simulate the next occurrence
### time using MCMC


poisfit <- function(faulti,inter.NA, t.occ, N, K,
                    n.iter.mc = 5010000,n.burnin.mc = 10000,n.thin.mc=1000){
	library(lattice)
	library(R2jags)	
	jagsdata <- list("inter.NA", "t.occ", "N", "K")
	
	# Define the parameters that we are interested in. 
	# Their posterior distributions will be the output.
	# t.occ.N is the forecast time of the next earthquake occurrence
	params <- c("lambda","sigmaZ","t.occ.N","tfore.mean","Z")

	# Set the initial values
	lam.est <- 1/mean(inter.NA[,-N])
	inits <- function(){
			list("lambda"=lam.est)
			}

	mod <- jags(data = jagsdata, inits = inits,
		            parameters.to.save = params, n.chains = 3,
	            n.iter = n.iter.mc, n.burnin = n.burnin.mc,n.thin=n.thin.mc,
	            model.file = textConnection(poisNoCensor.jags))
	
	# Convert to an MCMC object
	mod.mcmc <- as.mcmc(mod)
	
	### Can use the code below to save the result in a .image file		
	#	save(mod.mcmc,file=paste("../../Results/PoisSimRes2/PoisFault-mcmc",faulti,".image",sep=""))
	
	return(mod.mcmc)
	}
	
	



