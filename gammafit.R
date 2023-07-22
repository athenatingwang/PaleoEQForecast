# Fitting a gamma renewal process to each MC sample of the 
# earthquake chronologies from the output of an OxCal model

# This file contains the general jags run for any data.
# Refer to gammafit.allfaults.R which sources this gammafit.R and fits
# the Poisson process to the occurrence times from each fault.


tmp <- "
model{
  for (j in 1:K){ # jth MC sample of the data, K samples in total
   for (i in 1:N){ # ith observation of the data	
     # Likelihood 
     isCensor[j,i] ~ dinterval( inter.NA[j,i] , CensLim[j,i] )
     inter.NA[j,i] ~ dgamma(Z[j]*alpha,Y[j]*beta)  
     # Nested indexing for multiple samples
   }

   Z[j] ~ dgamma(sh, rt)
   Y[j] ~ dgamma(shy, rty)

   t.occ.N[j] <- t.occ[j,N] + inter.NA[j,N]  
  }

  beta ~ dnorm(0,10^-4)T(0,)
  alpha ~ dnorm(0,10^-4)T(0,)
  sh <- 1/pow(sigmaZ,2)
  rt <- 1/pow(sigmaZ,2)
  sigmaZ ~ dt(0,0.04,3)T(0,)
#  sigmaZ ~ dgamma(10,10)
  shy <- 1/pow(sigmaY,2)
  rty <- 1/pow(sigmaY,2)
  sigmaY ~ dt(0,0.04,3)T(0,)
#  sigmaY ~ dgamma(10,10)

  #  sigmaZ ~ dgamma(10,10)
  # t.occ.N is the forecast time of the next earthquake occurrence
  tfore.mean <- mean(t.occ.N)
#  tfore.med <- med(t.occ.N)
}
"

#curve(dgamma(x,shape=10,rate=10),0,10)

cat(tmp, file = "gamma.jags")

# faulti: integer indicating the ith fault in the folder
# chronologies_all
gammafit <- function(faulti,inter.NA, t.occ, N, K, isCensor, CensLim, plot.Y=FALSE){
	library(lattice)
	library(R2jags)	
	jagsdata <- list("inter.NA", "t.occ", "N", "K","isCensor","CensLim")
	
	# Define the parameters that we are interested in. 
	# Their posterior distributions will be the output.
	# t.occ.N is the forecast time of the next earthquake occurrence
	if (plot.Y){
	   params <- c("alpha","beta","sigmaY","sigmaZ","t.occ.N","tfore.mean","Y","Z")
	   }else{
	   params <- c("alpha","beta","sigmaY","sigmaZ","tfore.mean","Y","Z")
	   }

	# Set the initial values
	rate0 <- mean(inter.NA[,-N])/(mean(inter.NA[,-N]^2)-
	                       mean(inter.NA[,-N])^2)
	shape0 <- rate0*mean(inter.NA[,-N])
	
	inits <- function(){
			list("alpha"=shape0+runif(1,-shape0/4,shape0/4),"beta"=rate0+runif(1,-rate0/3,rate0/3))
			}

#	inter1 <- CensLim + 50
#	inter1[,-N] <- NA
#	inter2 <- inter1+50
#	inter3 <- inter2+50
	
#	inits <- function(){
#			list( list("alpha"=shape0,"beta"=rate0),
#			 #"inter.NA[,N]"=inter1[,N]),
#			 list("alpha"=shape0+runif(1,0,0.03),"beta"=rate0+runif(1,0,0.03)),
#			 #  "inter.NA[,N]"=inter2[,N]),
#			 list("alpha"=shape0+runif(1,0,0.03),"beta"=rate0+runif(1,0,0.03)) )
#			 #  "inter.NA[,N]"=inter3[,N]))
#			}


	mod <- jags(data = jagsdata, inits = inits,
		            parameters.to.save = params, n.chains = 3,
			    n.iter = 5010000, n.burnin = 10000,n.thin=1000,
			  #  n.iter = 1000, n.burnin = 100,n.thin=10,
			    model.file = 'gamma.jags')
	
	# Convert to an MCMC object
	mod.mcmc <- as.mcmc(mod)
	
	save(mod,file=paste("/nesi/nobackup/nesi00326/GammaRes/GammaFault",faulti,".image",sep=""))
	save(mod.mcmc,file=paste("/nesi/nobackup/nesi00326/GammaRes/GammaFault-mcmc",faulti,".image",sep=""))

	#plot(mod.mcmc[,c(1:6,6+K)],trace=TRUE)
	
	return(mod.mcmc)
	}
	
	







