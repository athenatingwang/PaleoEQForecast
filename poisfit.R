# Fitting a Poisson process to each MC sample of the earthquake chronologies
# from the output of an OxCal model

# This file contains the general jags run for any data.
# Refer to poisfit.allfaults.R which sources this poisfit.R and fits
# the Poisson process to the occurrence times from each fault.


tmp <- "
model{
  for (j in 1:K){ # jth MC sample of the data, K samples in total
   for (i in 1:N){ # ith observation of the data	
     # Likelihood 
     isCensor[j,i] ~ dinterval( inter.NA[j,i] , CensLim[j,i] )
     inter.NA[j,i] ~ dexp(lambda*Z[j]) # Nested indexing for multiple samples
   }
   Z[j] ~ dgamma(sh, rt)	
   t.occ.N[j] <- t.occ[j,N] + inter.NA[j,N]  
}
  sh <- 1/pow(sigmaZ,2)
  rt <- 1/pow(sigmaZ,2)
  sigmaZ ~ dt(0,0.04,3)T(0,)
#  sigmaZ ~ dgamma(10,10)
#  sigmaZ ~ dunif(0,10)
  lambda ~ dnorm(0,10^-4)T(0,)
  # t.occ.N is the forecast time of the next earthquake occurrence
  tfore.mean <- mean(t.occ.N)
}
"

#curve(dgamma(x,shape=10,rate=10),0,10)

cat(tmp, file = "pois.jags")

# faulti: integer indicating the ith fault in the folder
# chronologies_all
poisfit <- function(faulti,inter.NA, t.occ, N, K, isCensor, CensLim, plot.Y=TRUE){
	library(lattice)
	library(R2jags)	
	jagsdata <- list("inter.NA", "t.occ", "N", "K","isCensor","CensLim")
	
	# Define the parameters that we are interested in. 
	# Their posterior distributions will be the output.
	# t.occ.N is the forecast time of the next earthquake occurrence
	if (plot.Y){
	   params <- c("lambda","sigmaZ","t.occ.N","tfore.mean","Z")
	   }else{
	   params <- c("lambda","sigmaZ","tfore.mean")
	   }

	# Set the initial values
	lam.est <- 1/mean(inter.NA[,-N])
	inits <- function(){
			list("lambda"=lam.est+runif(1,-lam.est/4,lam.est/4))
			}

	mod <- jags(data = jagsdata, inits = inits,
		            parameters.to.save = params, n.chains = 3,
			    n.iter = 5010000, n.burnin = 10000, n.thin=1000,
			    #n.iter = 1000, n.burnin = 100, n.thin=5,
                            model.file = 'pois.jags')
	
	# Convert to an MCMC object
	mod.mcmc <- as.mcmc(mod)
	
	save(mod,file=paste("/nesi/nobackup/nesi00326/PoisRes/PoisFault",faulti,".image",sep=""))
	save(mod.mcmc,file=paste("/nesi/nobackup/nesi00326/PoisRes/PoisFault-mcmc",faulti,".image",sep=""))

       # if (plot.Y)
       #   plot(mod.mcmc[,c(1:4,4+K)],trace=TRUE)
	
	return(mod.mcmc)
	}
	
	














