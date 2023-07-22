# Fitting a lognormal renewal process to each MC sample of the 
# earthquake chronologies from the output of an OxCal model

# This file contains the general jags run for any data.
# Refer to lnormfit.allfaults.R which sources this lnormfit.R and fits
# the lognormal renewal process to the occurrence times from each fault.


tmp <- "
model{
  for (j in 1:K){ # jth MC sample of the data, K samples in total
   for (i in 1:N){ # ith observation of the data	
     # Likelihood 
     isCensor[j,i] ~ dinterval( inter.NA[j,i] , CensLim[j,i] )
     inter.NA[j,i] ~ dlnorm(Z[j]+mu,1/(Y[j]*sig)/(Y[j]*sig))  
     # Nested indexing for multiple samples
   }

   Z[j] ~ dnorm(0, 1/sigmaZ/sigmaZ)
   Y[j] ~ dgamma(shy, rty)

   t.occ.N[j] <- t.occ[j,N] + inter.NA[j,N]  
  }

  mu ~ dnorm(0,10^-4)
  sig ~ dt(0,0.04,3)T(0,)
  sigmaZ ~ dt(0,0.04,3)T(0,)
#  sigmaZ ~ dgamma(1,1)
  shy <- 1/pow(sigmaY,2)
  rty <- 1/pow(sigmaY,2)
  sigmaY ~ dt(0,0.04,3)T(0,)
#  sigmaY ~ dgamma(1,1)

  #  sigmaZ ~ dgamma(10,10)
  # t.occ.N is the forecast time of the next earthquake occurrence
  tfore.mean <- mean(t.occ.N)
#  tfore.med <- med(t.occ.N)
}
"

#curve(dgamma(x,shape=10,rate=10),0,10)

cat(tmp, file = "lnorm.jags")

# faulti: integer indicating the ith fault in the folder
# chronologies_all
lnormfit <- function(faulti,inter.NA, t.occ, N, K, isCensor, CensLim, plot.Y=TRUE){
	library(lattice)
	library(R2jags)	
	jagsdata <- list("inter.NA", "t.occ", "N", "K","isCensor","CensLim")
	
	# Define the parameters that we are interested in. 
	# Their posterior distributions will be the output.
	# t.occ.N is the forecast time of the next earthquake occurrence
	if (plot.Y){
	   params <- c("mu","sig","sigmaY","sigmaZ","t.occ.N","tfore.mean","Y","Z")
	   }else{
	   params <- c("mu","sig","sigmaY","sigmaZ","tfore.mean","Y","Z")
	   }

	# Set the initial values
	mu0 <- mean(log(inter.NA[,-N]))
	sig0 <- sd(log(inter.NA[,-N]))
	
	inits <- function(){
			list("mu"=mu0+runif(1,-mu0/5,mu0/5),"sig"=sig0+runif(1,-sig0/5,sig0/5))
			}

#	inits <- function(){#list("mu"=mu0,"sig"=sig0)
#			list( list("mu"=mu0,"sig"=sig0),
#			# "inter.NA[,N]"=inter1[,N]),
#			 list("mu"=mu0+runif(1,0.01,0.6),"sig"=sig0+runif(1,0.01,0.06)),
#			#   "inter.NA[,N]"=inter2[,N]),
#			 list("mu"=mu0+runif(1,0.01,0.6),"sig"=sig0+runif(1,0.01,0.06)) )
#			#   "inter.NA[,N]"=inter3[,N]))
#			}


	mod <- jags(data = jagsdata, inits = inits,
		            parameters.to.save = params, n.chains = 3,
			    n.iter = 5010000, n.burnin = 10000,n.thin=1000,
			   # n.iter = 1000, n.burnin = 100,n.thin=10,
			    model.file = 'lnorm.jags')
	
	# Convert to an MCMC object
	mod.mcmc <- as.mcmc(mod)
	

#	plot(mod.mcmc[,c(1:7,7+K)],trace=TRUE)
	
	return(list(mod=mod,mod.mcmc=mod.mcmc))
	}
	
	







