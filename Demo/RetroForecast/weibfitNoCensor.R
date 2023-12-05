# Fitting a Weibull renewal process without censoring to each MC sample of the earthquake chronologies 

# This file contains the jags model and general jags run for any data.
# Refer to weibfit.allfaults.R which sources this weibfitNoCensor.R and fits
# the Weibull renewal process without censoring to the occurrence times from each fault segment.


###############################################################
############               Jags model              ############
###############################################################
weibNoCensor.jags <- "
model{
  for (j in 1:K){ # jth MC sample of the data, K samples in total
   for (i in 1:N){ # ith observation of the jth MC sample	
     # Likelihood 
     inter.NA[j,i] ~ dweib(Z[j]*alpha,pow(Y[j]*beta,Z[j]*alpha))  
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
#### Refer to weibfit.allfaults.R for how to get each input argument below from the 
#### earthquake chronologies.
###################################################################################
# faulti: integer indicating the ith fault in the folder chronologies_all_final
# N: Number of earthquakes at the fault
# K: Number of MC samples used from the earthquake chronologies (we used 100)

# inter.NA: KxN matrix. It contains the interevent times from the K MC samples, 
### and the last column is NA as the next earthquake hasn't occurred and we
### want to forecast its occurrence time. JAGS will simulate the next occurrence
### time using MCMC


weibfit <- function(faulti,inter.NA, t.occ, N, K,
                    n.iter.mc = 5010000,n.burnin.mc = 10000,n.thin.mc=1000){
	library(lattice)
	library(R2jags)	
	jagsdata <- list("inter.NA", "t.occ", "N", "K")
	
	# Define the parameters that we are interested in. 
	# Their posterior distributions will be the output.
	# t.occ.N is the forecast time of the next earthquake occurrence
	params <- c("alpha","beta","sigmaY","sigmaZ","t.occ.N","tfore.mean","Y","Z")

	# Set the initial values for Weibull
	interarrival <- inter.NA[1:5,-N]
	n <- length(interarrival)
	yi <- zi <- NULL
	for (jj in 1:n){
	   yi[jj] <- log( quantile(interarrival,jj/(n+1)) )
	   zi[jj] <- log( -log(1-jj/(n+1)) )
	}
	medb1i <- medb0i <- NULL
	for (k in 1:n){
	   medb1i[k] <- median((yi[k]-yi[-k])/(zi[k]-zi[-k]))
	   medb0i[k] <- median((zi[k]*yi[-k]-zi[-k]*yi[k])/(zi[k]-zi[-k]))
	}
	b1 <- median(medb1i)
	if(b1==0) b1 <- mean(medb1i)
	b0 <- median(medb0i)
	shape0 <- 1/b1
	rate0 <- 1/exp(b0)

	inits <- function(){
			list("alpha"=shape0+runif(1,-shape0/5,shape0/5),"beta"=rate0+runif(1,-rate0/5,rate0/5))
			}


	mod <- jags(data = jagsdata, inits = inits,
		            parameters.to.save = params, n.chains = 3,
	            n.iter = n.iter.mc, n.burnin = n.burnin.mc,n.thin=n.thin.mc,
			    model.file = textConnection(weibNoCensor.jags))
	
	# Convert to an MCMC object
	mod.mcmc <- as.mcmc(mod)
	
	### Can use the code below to save the result in a .image file		
	save(mod.mcmc,file=paste("../../ResultsDemo/WeibSimRes2/WeibFault-mcmc",faulti,".image",sep=""))
	
	return(mod.mcmc)
	}
	


	