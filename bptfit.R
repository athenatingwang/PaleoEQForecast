# Fitting a BPT renewal process to each MC sample of the 
# earthquake chronologies from the output of an OxCal model

# This file contains the general jags run for any data.
# Refer to bptfit.allfaults.R which sources this bptfit.R and fits
# the BPT renewal process to the occurrence times from each fault.


tmp <- "
# Use parameterisation from Matthews (2002) BSSA paper
data {
 for (j in 1:K){
   for (i in 1:N){
     zeros[j,i] = 0
   }
 }
 C = 10000
 PI = 3.14159265358979323846
}

model{
  for (j in 1:K){ # jth MC sample of the data, K samples in total
    # Censored times
    # survival function for the censored time
    #isCensor[j,N] ~ dinterval( inter.last[j] , CensLim[j,N] )
    #isCensor[j,N] ~ dbern(1-cens.prob[j])
    #inter.last[j] ~ dunif(CensLim[j,N]+0.5,10000000)
    zeros[j,N] ~ dpois(zeros.mean[j,N])
    zeros.mean[j,N] = -l[j,N] + C

    shape[j] = alpha * Y[j]
    loc[j] = mu * Z[j]

    #l[j,N] = 0.5 * (log(loc[j]) - log(2*PI) - 2*log(shape[j]) - 
    #   3*log(inter.last[j])) - pow((inter.last[j] - loc[j]), 2)/
    #   (2*loc[j]*pow(shape[j],2)*inter.last[j])

    cens1[j] = (1/shape[j]) * pow(loc[j]/CensLim[j,N], 1/2) *
    	(CensLim[j,N]/loc[j] - 1)
    cens2[j] = -(1/shape[j]) * pow(loc[j]/CensLim[j,N], 1/2) *
    	(CensLim[j,N]/loc[j] + 1)
    cens.prob[j] = pnorm(cens1[j],0,1) + 
      exp(2/(pow(shape[j], 2)))*pnorm(cens2[j],0,1)

    l[j,N] = log(1-cens.prob[j])
	
    ## likelihood using zero trick
    for (i in 1:(N-1)){ # ith observation of the data
      #isCensor[j,i] ~ dinterval( inter.NA[j,i] , CensLim[j,i] )
      zeros[j,i] ~ dpois(zeros.mean[j,i])
      zeros.mean[j,i] = -l[j,i] + C
      # Log-likelihood
      l[j,i] = 0.5 * (log(loc[j]) - log(2*PI) - 2*log(shape[j]) - 3*log(inter.NA[j,i])) -
       pow((inter.NA[j,i] - loc[j]), 2)/(2*loc[j]*pow(shape[j],2)*inter.NA[j,i])
    }
  	Z[j] ~ dgamma(shz, rtz)
    Y[j] ~ dgamma(shy, rty)

    #t.occ.N[j] <- t.occ[j,N] + inter.last[j]  
  }

  mu ~ dnorm(0,10^-4)T(0,)
  alpha ~ dt(0,0.04,3)T(0,)
  shz <- 1/pow(sigmaZ,2)
  rtz <- 1/pow(sigmaZ,2)
  sigmaZ ~ dt(0,0.04,3)T(0,)
  shy <- 1/pow(sigmaY,2)
  rty <- 1/pow(sigmaY,2)
  sigmaY ~ dt(0,0.04,3)T(0,)
#  sigmaY ~ dgamma(10,10)
  # t.occ.N is the forecast time of the next earthquake occurrence
  #tfore.mean <- mean(t.occ.N)
}
"

#curve(dgamma(x,shape=10,rate=10),0,10)

cat(tmp, file = "bpt.jags")

# faulti: integer indicating the ith fault in the folder
# chronologies_all
#bptfit <- function(faulti,inter.NA, t.occ, N, K, isCensor, CensLim, plot.Y=FALSE){
bptfit <- function(faulti,inter.NA, N, K, CensLim, plot.Y=FALSE){
  library(lattice)
	library(R2jags)	
	#library(actuar)
	#jagsdata <- list("inter.NA", "t.occ", "N", "K","isCensor","CensLim")
  jagsdata <- list("inter.NA", "N", "K", "CensLim")
  
	# Define the parameters that we are interested in. 
	# Their posterior distributions will be the output.
	# t.occ.N is the forecast time of the next earthquake occurrence
	if (plot.Y){
	   params <- c("mu","alpha","sigmaY","sigmaZ","Y","Z")#,"t.occ.N","tfore.mean")#,"tfore.med")
	   }else{
	   params <- c("mu","alpha","sigmaY","sigmaZ","Y","Z")#,"tfore.mean")#,"tfore.med")
	   }

	# Set the initial values
	mu0 <- mean(inter.NA[,-N])
	alpha0 <- max(sqrt(mean(1/inter.NA[,-N])*mu0-1),0.0001)
	

#	inter1 <- CensLim + 50
#	inter1[,-N] <- NA
#	inter2 <- inter1+50
#	inter3 <- inter2+50
	
	inits <- function(){list("mu"=mu0+runif(1,-mu0/4,mu0/4),"alpha"=alpha0+runif(1,-alpha0/3,alpha0/3))
#			list( list("mu"=mu0,"alpha"=alpha0),
			# "inter.NA[,N]"=inter1[,N]),
#			 list("mu"=mu0+runif(1,0.01,0.06),"alpha"=alpha0+runif(1,0.01,0.06)),
			#   "inter.NA[,N]"=inter2[,N]),
#			 list("mu"=mu0+runif(1,0.01,0.06),"alpha"=alpha0+runif(1,0.01,0.06)) )
			#   "inter.NA[,N]"=inter3[,N]))
			}

	mod <- jags(data = jagsdata, inits = inits,
		            parameters.to.save = params, n.chains = 3,
			    n.iter = 5010000, n.burnin = 10000,n.thin=1000,
			  #  n.iter = 1000, n.burnin = 100,n.thin=10,
			    model.file = 'bpt.jags')
	
	# Convert to an MCMC object
	mod.mcmc <- as.mcmc(mod)
	
	save(mod,file=paste("/nesi/nobackup/nesi00326/bptRes/bptFault",faulti,".image",sep=""))
	save(mod.mcmc,file=paste("/nesi/nobackup/nesi00326/bptRes/bptFault-mcmc",faulti,".image",sep=""))

#	plot(mod.mcmc[,c(1:7,7+K)],trace=TRUE)
	
	return(mod.mcmc)
	}
	
	







