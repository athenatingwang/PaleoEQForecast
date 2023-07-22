# Fitting a Weibull renewal process to each MC sample of the 
# earthquake chronologies from the output of an OxCal model

# This file contains the general jags run for any data.
# Refer to weibfit.allfaults.R which sources this weibfit.R and fits
# the Weibull renewal process to the occurrence times from each fault.


tmp <- "
model{
  for (j in 1:K){ # jth MC sample of the data, K samples in total
   for (i in 1:N){ # ith observation of the data	
     # Likelihood 
     isCensor[j,i] ~ dinterval( inter.NA[j,i] , CensLim[j,i] )
     inter.NA[j,i] ~ dweib(Z[j]*alpha,pow(Y[j]*beta,Z[j]*alpha))  
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

cat(tmp, file = "weib.jags")

# faulti: integer indicating the ith fault in the folder
# chronologies_all
weibfit <- function(faulti,inter.NA, t.occ, N, K, isCensor, CensLim, plot.Y=FALSE){
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
			list("alpha"=shape0+runif(1,-shape0/3,shape0/3),"beta"=rate0+runif(1,-rate0/3,rate0/3))
			}

#	inits <- function(){
#			list( list("alpha"=shape0,"beta"=1/rate0),
#			 #"inter.NA[,N]"=inter1[,N]),
#			 list("alpha"=shape0+runif(1,0,0.03),"beta"=1/rate0+runif(1,0,0.03)),
#			 #  "inter.NA[,N]"=inter2[,N]),
#			 list("alpha"=shape0+runif(1,0,0.03),"beta"=1/rate0+runif(1,0,0.03)) )
#			 #  "inter.NA[,N]"=inter3[,N]))
#			}

	mod <- jags(data = jagsdata, inits = inits,
		            parameters.to.save = params, n.chains = 3,
			    n.iter = 5010000, n.burnin = 10000,n.thin=1000,
			    model.file = 'weib.jags')
	
	# Convert to an MCMC object
	mod.mcmc <- as.mcmc(mod)
	
	save(mod,file=paste("/nesi/nobackup/nesi00326/WeibRes/WeibFault",faulti,".image",sep=""))
	save(mod.mcmc,file=paste("/nesi/nobackup/nesi00326/WeibRes/WeibFault-mcmc",faulti,".image",sep=""))


	#plot(mod.mcmc[,c(1:6,6+K)],trace=TRUE)
	
	return(mod.mcmc)
	}
	
	