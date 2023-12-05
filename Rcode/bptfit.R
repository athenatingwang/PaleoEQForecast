# Fitting a BPT renewal process to each MC sample of the earthquake chronologies 

# This file contains the jags model and general jags run for any data.
# Refer to bptfit.allfaults.R which sources this bptfit.R and fits
# the BPT renewal process to the occurrence times from each fault segment.


###############################################################
############               Jags model              ############
###############################################################
bpt.jags <- "
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
    zeros[j,N] ~ dpois(zeros.mean[j,N])
    zeros.mean[j,N] = -l[j,N] + C

    shape[j] = alpha * Y[j]
    loc[j] = mu * Z[j]

    cens1[j] = (1/shape[j]) * pow(loc[j]/CensLim[j,N], 1/2) *
    	(CensLim[j,N]/loc[j] - 1)
    cens2[j] = -(1/shape[j]) * pow(loc[j]/CensLim[j,N], 1/2) *
    	(CensLim[j,N]/loc[j] + 1)
    cens.prob[j] = pnorm(cens1[j],0,1) + 
      exp(2/(pow(shape[j], 2)))*pnorm(cens2[j],0,1)

    l[j,N] = log(1-cens.prob[j])
	
    for (i in 1:(N-1)){ # ith observation of the jth MC sample
      zeros[j,i] ~ dpois(zeros.mean[j,i])
      zeros.mean[j,i] = -l[j,i] + C
      # Log-likelihood
      l[j,i] = 0.5 * (log(loc[j]) - log(2*PI) - 2*log(shape[j]) - 3*log(inter.NA[j,i])) -
       pow((inter.NA[j,i] - loc[j]), 2)/(2*loc[j]*pow(shape[j],2)*inter.NA[j,i])
    }
    Z[j] ~ dgamma(shz, rtz)
    Y[j] ~ dgamma(shy, rty)
  }

  mu ~ dnorm(0,10^-4)T(0,)
  alpha ~ dt(0,0.04,3)T(0,)
  shz <- 1/pow(sigmaZ,2)
  rtz <- 1/pow(sigmaZ,2)
  sigmaZ ~ dt(0,0.04,3)T(0,)
  shy <- 1/pow(sigmaY,2)
  rty <- 1/pow(sigmaY,2)
  sigmaY ~ dt(0,0.04,3)T(0,)
}
"
###############################################################
############               Jags model              ############
###############################################################




###################################################################################
#### Refer to bptfit.allfaults.R for how to get each input argument below from the 
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


bptfit <- function(faulti,inter.NA, N, K, CensLim,
                   n.iter.mc = 5010000,n.burnin.mc = 10000,n.thin.mc=1000){
  library(lattice)
  library(R2jags)	
  jagsdata <- list("inter.NA", "N", "K", "CensLim")
  
  # Define the parameters that we are interested in. 
  # Their posterior distributions will be the output.
  ### For BPT renewal process, t.occ.N, the forecast time of the next earthquake occurrence
  ### is obtained after fitting the model using the posterior samples of the parameters.
  ### The forecast code is in bptfit.allfaults.R
  params <- c("mu","alpha","sigmaY","sigmaZ","Y","Z")

  # Set the initial values
  mu0 <- mean(inter.NA[,-N])
  alpha0 <- max(sqrt(mean(1/inter.NA[,-N])*mu0-1),0.0001)
	
	
  inits <- function(){
             list("mu"=mu0+runif(1,-mu0/4,mu0/4),"alpha"=alpha0+runif(1,-alpha0/3,alpha0/3))
		}

  mod <- jags(data = jagsdata, inits = inits,
		  parameters.to.save = params, n.chains = 3,
	        n.iter = n.iter.mc, n.burnin = n.burnin.mc,n.thin=n.thin.mc,
	        model.file = textConnection(bpt.jags))
	
  # Convert to an MCMC object
  mod.mcmc <- as.mcmc(mod)
	
### Can use the code below to save the result in a .image file
  save(mod.mcmc,file=paste("../Results/bptRes/bptFault-mcmc",faulti,".image",sep=""))
	
  return(mod.mcmc)
}
	
	







