library(rjags)
a <- dir("../DataFinal/chronologies_all_final")
nfault <- length(a)

K = 100

Noccur <- NULL
waic <- weights <- matrix(NA,nrow=nfault,ncol=5)

p50Pois <- p100Pois <- p500Pois <- array(NA,c(nfault,K,15000))
p50gam <- p100gam <- p500gam <- array(NA,c(nfault,K,15000))
p50weib <- p100weib <- p500weib <- array(NA,c(nfault,K,15000))
p50bpt <- p100bpt <- p500bpt <- array(NA,c(nfault,K,15000))
p50lnorm <- p100lnorm <- p500lnorm <- array(NA,c(nfault,K,15000))

p50modAve <- p100modAve <- p500modAve <- array(NA,c(nfault,K,15000))

modAveFore <- Fore.Pois <- Fore.gam <- matrix(NA,nrow=nfault,ncol=1500000)
Fore.weib <- Fore.bpt <- Fore.lnorm <- matrix(NA,nrow=nfault,ncol=1500000)

for (faulti in (1:nfault)){
  # Load data
  data = read.csv(paste('../DataFinal/chronologies_all_final/',a[faulti],sep=''), header=FALSE)
  
  t.occ = data.matrix(data)
  m = nrow(t.occ)
  N <- ncol(t.occ)
  Noccur[faulti] <- N  
  # Add a column with 2022 to indicate
  # the censored time at the year 2022, for forecasting purpose
  t.occ.cens = cbind(t.occ,rep(2022,m))
  # Calculate inter-event times. 
  # The last column is the censored time 
  inter.cens = t(diff(t(t.occ.cens))) # Transpose, take difference and transpose back


  # Poisson process results
  load(paste("../Results/PoisRes/PoisFault-mcmc",faulti,".image",sep=""))

  pois.mat <- as.matrix(mod.mcmc)
  pois.dat <- as.data.frame(pois.mat)
  lambda <- pois.dat$lambda
  n.mc <- length(lambda)
  Fore.Pois[faulti,] <- as.vector(pois.mat[,4:(3+K)])

  
  pPois <- plogPois <- array(NA,c(K,N,length(lambda)))
  ppdPois <- array(NA,c(K,N))
  pwaicPois <- array(NA,c(K,N))

  hpois <- function(x,rate){
    dexp(x,rate)/pexp(x,rate,lower.tail = FALSE)
  }   
  Fcondpois <- function(rate,lower,upper){
    1-exp(-integrate(hpois,lower=lower,upper=upper,rate=rate)$value)
  }
  
  for (j in 1:K){ # jth MC sample of the data, K samples in total
    Zj <- eval(parse(text=paste('pois.dat$"Z[',j,']"',sep="")))
    for (i in 1:(N-1)){ # ith observation of the data	
      pPois[j,i,] <- dexp(inter.cens[j,i],lambda*Zj)
      ppdPois[j,i] <- mean(pPois[j,i,])
      plogPois[j,i,] <- dexp(inter.cens[j,i],lambda*Zj,log=TRUE)
      pwaicPois[j,i] <- var(plogPois[j,i,])
    }
    pPois[j,N,] <- pexp(inter.cens[j,N],lambda*Zj, lower.tail = FALSE)
    ppdPois[j,N] <- mean(pPois[j,N,])
    plogPois[j,N,] <- pexp(inter.cens[j,N],lambda*Zj,lower.tail=FALSE,log.p=TRUE)
    pwaicPois[j,N] <- var(plogPois[j,N,])

    x50 <- inter.cens[j,N]+50
    x100 <- inter.cens[j,N]+100
    x500 <- inter.cens[j,N]+500
    x.cens <- inter.cens[j,N]

	p50Pois[faulti,j,] <- apply(t(lambda*Zj),2,Fcondpois,lower=x.cens,upper=x50)
    p100Pois[faulti,j,] <- apply(t(lambda*Zj),2,Fcondpois,lower=x.cens,upper=x100)
    p500Pois[faulti,j,] <- apply(t(lambda*Zj),2,Fcondpois,lower=x.cens,upper=x500)
}

  sum(log(ppdPois))
  sum(pwaicPois)

  waicPois <- -2*sum(log(ppdPois)) + 2* sum(pwaicPois)



  # Gamma renewal process results
  load(paste("../Results/GammaRes/GammaFault-mcmc",faulti,".image",sep=""))

  gam.mat <- as.matrix(mod.mcmc)
  gam.dat <- as.data.frame(gam.mat)
  alpha <- gam.dat$alpha
  beta <- gam.dat$beta
  n.mc <- length(alpha)
  Fore.gam[faulti,] <- as.vector(gam.mat[,6:(5+K)])
  
  
  pgam <- ploggam <- array(NA,c(K,N,length(alpha)))
  ppdgam <- array(NA,c(K,N))
  pwaicgam <- array(NA,c(K,N))
 
  hgam <- function(x,param){
    shape <- param[1]
    rate <- param[2]
    dgamma(x,shape,rate)/pgamma(x,shape,rate,lower.tail = FALSE)
  }   
  Fcondgam <- function(param,lower,upper){
    1-exp(-integrate(hgam,lower=lower,upper=upper,param=param)$value)
  }
 
  for (j in 1:K){ # jth MC sample of the data, K samples in total
    Yj <- eval(parse(text=paste('gam.dat$"Y[',j,']"',sep="")))
    Zj <- eval(parse(text=paste('gam.dat$"Z[',j,']"',sep="")))
    for (i in 1:(N-1)){ # ith observation of the data	
      pgam[j,i,] <- dgamma(inter.cens[j,i],shape=Zj*alpha,rate=Yj*beta)
      ppdgam[j,i] <- mean(pgam[j,i,])
      ploggam[j,i,] <- dgamma(inter.cens[j,i],shape=Zj*alpha,rate=Yj*beta,log=TRUE)
      pwaicgam[j,i] <- var(ploggam[j,i,])
    }
    pgam[j,N,] <- pgamma(inter.cens[j,N],shape=Zj*alpha,rate=Yj*beta,lower.tail=FALSE)
    ppdgam[j,N] <- mean(pgam[j,N,])
    ploggam[j,N,] <- pgamma(inter.cens[j,N],shape=Zj*alpha,rate=Yj*beta,lower.tail=FALSE,log.p=TRUE)
    pwaicgam[j,N] <- var(ploggam[j,N,])

    x50 <- inter.cens[j,N]+50
    x100 <- inter.cens[j,N]+100
    x500 <- inter.cens[j,N]+500
    x.cens <- inter.cens[j,N]

    params <- cbind(Zj*alpha,Yj*beta)
	for (i in 1:15000){
      p50gam[faulti,j,i] <- Fcondgam(params[i,],lower=x.cens,upper=x50)
      p100gam[faulti,j,i] <- Fcondgam(params[i,],lower=x.cens,upper=x100)
      p500gam[faulti,j,i] <- Fcondgam(params[i,],lower=x.cens,upper=x500)
	}
}
  
  sum(log(ppdgam))
  sum(pwaicgam)
  
  waicgam <- -2*sum(log(ppdgam)) + 2* sum(pwaicgam)
  
  
  
  # Weibull renewal process results
  load(paste("../Results/WeibRes/WeibFault-mcmc",faulti,".image",sep=""))

  weib.mat <- as.matrix(mod.mcmc)
  weib.dat <- as.data.frame(weib.mat)
  alpha <- weib.dat$alpha
  beta <- weib.dat$beta
  n.mc <- length(alpha)
  Fore.weib[faulti,] <- as.vector(weib.mat[,6:(5+K)])
  
  
  pweib <- plogweib <- array(NA,c(K,N,length(alpha)))
  ppdweib <- array(NA,c(K,N))
  pwaicweib <- array(NA,c(K,N))
  
  hweib <- function(x,param){
    shape <- param[1]
    scale <- param[2]
    dweibull(x,shape,scale)/pweibull(x,shape,scale,lower.tail = FALSE)
  }   
  Fcondweib <- function(param,lower,upper){
    1-exp(-integrate(hweib,lower=lower,upper=upper,param=param)$value)
  }

  for (j in 1:K){ # jth MC sample of the data, K samples in total
    Yj <- eval(parse(text=paste('weib.dat$"Y[',j,']"',sep="")))
    Zj <- eval(parse(text=paste('weib.dat$"Z[',j,']"',sep="")))
    shapej <- Zj*alpha
    scalej <- (Yj*beta)^(-1)
    for (i in 1:(N-1)){ # ith observation of the data
      pweib[j,i,] <- dweibull(inter.cens[j,i],shape=shapej,scale=scalej)
      ppdweib[j,i] <- mean(pweib[j,i,])
      plogweib[j,i,] <- dweibull(inter.cens[j,i],shape=shapej,scale=scalej,log=TRUE)
      pwaicweib[j,i] <- var(plogweib[j,i,])
    }
    pweib[j,N,] <- pweibull(inter.cens[j,N],shape=shapej,scale=scalej,lower.tail=FALSE)
    ppdweib[j,N] <- mean(pweib[j,N,])
    plogweib[j,N,] <- pweibull(inter.cens[j,N],shape=shapej,scale=scalej,
                               lower.tail=FALSE,log.p=TRUE)
    pwaicweib[j,N] <- var(plogweib[j,N,])
    
    x50 <- inter.cens[j,N]+50
    x100 <- inter.cens[j,N]+100
    x500 <- inter.cens[j,N]+500
    x.cens <- inter.cens[j,N]

    params <- cbind(shapej,scalej)
	for (i in 1:15000){
      p50weib[faulti,j,i] <- Fcondweib(params[i,],lower=x.cens,upper=x50)
      p100weib[faulti,j,i] <- Fcondweib(params[i,],lower=x.cens,upper=x100)
      p500weib[faulti,j,i] <- Fcondweib(params[i,],lower=x.cens,upper=x500)
	}
  }
  
  sum(log(ppdweib))
  sum(pwaicweib)
  
  waicweib <- -2*sum(log(ppdweib)) + 2* sum(pwaicweib)
  
  

  # BPT renewal process results
  load(paste("../Results/bptRes/bptFault-mcmc",faulti,".image",sep=""))
  load(paste("../Results/bptRes/HPD-toccBPT",faulti,".image",sep=""))
  library(actuar)
  
  bpt.mat <- as.matrix(mod.mcmc)
  bpt.dat <- as.data.frame(bpt.mat)
  alpha <- bpt.dat$alpha
  mu <- bpt.dat$mu
  n.mc <- length(alpha)
  Fore.bpt[faulti,] <- as.vector(bpt.res$t.occ.N)


  pbpt <- plogbpt <- array(NA,c(K,N,length(alpha)))
  ppdbpt <- array(NA,c(K,N))
  pwaicbpt <- array(NA,c(K,N))

  hbpt <- function(x,param){
    mean <- param[1]
	shape <- param[2]
    dinvgauss(x,mean,shape)/pinvgauss(x,mean,shape,lower.tail = FALSE)
  }   
  Fcondbpt <- function(param,lower,upper){
    1-exp(-integrate(hbpt,lower=lower,upper=upper,param=param)$value)
  }

  for (j in 1:K){ # jth MC sample of the data, K samples in total
    Yj <- eval(parse(text=paste('bpt.dat$"Y[',j,']"',sep="")))
    Zj <- eval(parse(text=paste('bpt.dat$"Z[',j,']"',sep="")))
    alphaj = alpha * Yj
    muj = mu * Zj
    shapej <- muj/alphaj^2
    
    for (i in 1:(N-1)){ # ith observation of the data
      pbpt[j,i,] <- dinvgauss(inter.cens[j,i], mean=muj, shape = shapej)
      ppdbpt[j,i] <- mean(pbpt[j,i,])
      plogbpt[j,i,] <- dinvgauss(inter.cens[j,i], mean=muj, shape = shapej,log=TRUE)
      pwaicbpt[j,i] <- var(plogbpt[j,i,])
    }
    pbpt[j,N,] <- pinvgauss(inter.cens[j,N], mean=muj, shape = shapej,lower.tail=FALSE)
    ppdbpt[j,N] <- mean(pbpt[j,N,])
    plogbpt[j,N,] <- pinvgauss(inter.cens[j,N], mean=muj, shape = shapej,
                               lower.tail=FALSE,log.p=TRUE)
    pwaicbpt[j,N] <- var(plogbpt[j,N,])
    
    x50 <- inter.cens[j,N]+50
    x100 <- inter.cens[j,N]+100
    x500 <- inter.cens[j,N]+500
    x.cens <- inter.cens[j,N]

    params <- cbind(muj,shapej)
	for (i in 1:15000){
      p50bpt[faulti,j,i] <- Fcondbpt(params[i,],lower=x.cens,upper=x50)
      p100bpt[faulti,j,i] <- Fcondbpt(params[i,],lower=x.cens,upper=x100)
      p500bpt[faulti,j,i] <- Fcondbpt(params[i,],lower=x.cens,upper=x500)
	}
  }
  
  sum(log(ppdbpt))
  sum(pwaicbpt)
  
  waicbpt <- -2*sum(log(ppdbpt)) + 2* sum(pwaicbpt)
  
  
  
  # lognormal renewal process results
  load(paste("../Results/lnormRes/lnormFault-mcmc",faulti,".image",sep=""))

  lnorm.mat <- as.matrix(mod.mcmc)
  lnorm.dat <- as.data.frame(lnorm.mat)
  mu <- lnorm.dat$mu
  sig <- lnorm.dat$sig
  n.mc <- length(mu)
  Fore.lnorm[faulti,] <- as.vector(lnorm.mat[,6:(5+K)])
  

  plnorm <- ploglnorm <- array(NA,c(K,N,length(alpha)))
  ppdlnorm <- array(NA,c(K,N))
  pwaiclnorm <- array(NA,c(K,N))
  
  hlnorm <- function(x,param){
    meanlog <- param[1]
	sdlog <- param[2]
    exp(dlnorm(x,meanlog, sdlog,log=TRUE)-plnorm(x,meanlog, sdlog,log.p=TRUE,lower.tail = FALSE))
  }   
  Fcondlnorm <- function(param,lower,upper){
    1-exp(-integrate(hlnorm,lower=lower,upper=upper,param=param)$value)
  }

  for (j in 1:K){ # jth MC sample of the data, K samples in total
    Yj <- eval(parse(text=paste('lnorm.dat$"Y[',j,']"',sep="")))
    Zj <- eval(parse(text=paste('lnorm.dat$"Z[',j,']"',sep="")))
  	meanlj <- Zj+mu
	  sdlj <- Yj*sig
    for (i in 1:(N-1)){ # ith observation of the data	
      plnorm[j,i,] <- dlnorm(inter.cens[j,i],meanlj,sdlj)
      ppdlnorm[j,i] <- mean(plnorm[j,i,])
      ploglnorm[j,i,] <- dlnorm(inter.cens[j,i],meanlj,sdlj,log=TRUE)
      pwaiclnorm[j,i] <- var(ploglnorm[j,i,])
    }
    plnorm[j,N,] <- plnorm(inter.cens[j,N],meanlj,sdlj,lower.tail=FALSE)
    ppdlnorm[j,N] <- mean(plnorm[j,N,])
    ploglnorm[j,N,] <- plnorm(inter.cens[j,N],meanlj,sdlj,
                               lower.tail=FALSE,log.p=TRUE)
    pwaiclnorm[j,N] <- var(ploglnorm[j,N,])

    x50 <- inter.cens[j,N]+50
    x100 <- inter.cens[j,N]+100
    x500 <- inter.cens[j,N]+500
    x.cens <- inter.cens[j,N]

    params <- cbind(meanlj,sdlj)
	for (i in 1:15000){
      p50lnorm[faulti,j,i] <- Fcondlnorm(params[i,],lower=x.cens,upper=x50)
      p100lnorm[faulti,j,i] <- Fcondlnorm(params[i,],lower=x.cens,upper=x100)
      p500lnorm[faulti,j,i] <- Fcondlnorm(params[i,],lower=x.cens,upper=x500)
	}
  }
  
  sum(log(ppdlnorm))
  sum(pwaiclnorm)
  
  waiclnorm <- -2*sum(log(ppdlnorm)) + 2* sum(pwaiclnorm)
  
  
  waic[faulti,] <- c(waicPois,waicgam,waicweib,waicbpt,waiclnorm)
  dwaic <- waic[faulti,] - min(waic[faulti,])
  temp.wt <- exp(-dwaic / 2)
  weights[faulti,] <- temp.wt / sum(temp.wt)

  
  set.seed(123)
  ind <- sample(1:5,length(alpha),replace=T,prob=weights[faulti,])
  modAveFore[faulti,] <- (ind==1)*Fore.Pois[faulti,] + (ind==2)*Fore.gam[faulti,] +
    (ind==3)*Fore.weib[faulti,] + (ind==4)*Fore.bpt[faulti,] +
    (ind==5)*Fore.lnorm[faulti,]

  p50modAve[faulti,,] <- (ind==1)*p50Pois[faulti,,] + (ind==2)*p50gam[faulti,,] +
    (ind==3)*p50weib[faulti,,] + (ind==4)*p50bpt[faulti,,] +
    (ind==5)*p50lnorm[faulti,,] 
  p100modAve[faulti,,] <- (ind==1)*p100Pois[faulti,,] + (ind==2)*p100gam[faulti,,] +
    (ind==3)*p100weib[faulti,,] + (ind==4)*p100bpt[faulti,,] +
    (ind==5)*p100lnorm[faulti,,]
  p500modAve[faulti,,] <- (ind==1)*p500Pois[faulti,,] + (ind==2)*p500gam[faulti,,] +
    (ind==3)*p500weib[faulti,,] + (ind==4)*p500bpt[faulti,,] +
    (ind==5)*p500lnorm[faulti,,] 

print(faulti)

}


#ave.probres.all <- list("p50Pois"=p50Pois,"p100Pois"=p100Pois,"p500Pois"=p500Pois,
#              "p50gam"=p50gam,"p100gam"=p100gam,"p500gam"=p500gam,
#              "p50weib"=p50weib,"p100weib"=p100weib,"p500weib"=p500weib,
#              "p50bpt"=p50bpt,"p100bpt"=p100bpt,"p500bpt"=p500bpt,
#              "p50lnorm"=p50lnorm,"p100lnorm"=p100lnorm,"p500lnorm"=p500lnorm,
#              "p50modAve"=p50modAve,"p100modAve"=p100modAve,"p500modAve"=p500modAve) 

#ave.logprobres.all <- list("lp50Pois"=log(p50Pois),"lp100Pois"=log(p100Pois),"lp500Pois"=log(p500Pois),
#              "lp50gam"=log(p50gam),"lp100gam"=log(p100gam),"lp500gam"=log(p500gam),
#              "lp50weib"=log(p50weib),"lp100weib"=log(p100weib),"lp500weib"=log(p500weib),
#              "lp50bpt"=log(p50bpt),"lp100bpt"=log(p100bpt),"lp500bpt"=log(p500bpt),
#              "lp50lnorm"=log(p50lnorm),"lp100lnorm"=log(p100lnorm),"lp500lnorm"=log(p500lnorm),
#              "lp50modAve"=log(p50modAve),"lp100modAve"=log(p100modAve),"lp500modAve"=log(p500modAve)) 


ave.foreres <- list("modAveFore"=modAveFore,"Fore.Pois"=Fore.Pois,
                    "Fore.gam"=Fore.gam,"Fore.weib"=Fore.weib,
                    "Fore.bpt"=Fore.bpt,"Fore.lnorm"=Fore.lnorm,
                    "weights"=weights,"waic"=waic) 

ave.probres <- list("p50Pois"=p50Pois,"p50gam"=p50gam,"p50weib"=p50weib,
              "p50bpt"=p50bpt,"p50lnorm"=p50lnorm,"p50modAve"=p50modAve) 

ave.logprobres <- list("lp50Pois"=log(p50Pois),"lp50gam"=log(p50gam),"lp50weib"=log(p50weib),
              "lp50bpt"=log(p50bpt),"lp50lnorm"=log(p50lnorm),"lp50modAve"=log(p50modAve)) 

save(ave.probres,file=paste("ModAveProbRes.image",sep=""))
save(ave.logprobres,file=paste("ModAveLogProbRes.image",sep=""))
save(ave.foreres,file=paste("ModAveForeRes.image",sep=""))

write.csv(ave.foreres$weights,"WAICweights.csv")








pois.lambda <- gam.alpha <- gam.beta <- weib.alpha <- weib.beta <- matrix(NA,nfault,3)
lnorm.mu <- lnorm.sig <- bpt.alpha <- bpt.mu <- matrix(NA,nfault,3)

for (faulti in 1:nfault){
  # Poisson process results
  load(paste("../Results/PoisRes/PoisFault-mcmc",faulti,".image",sep=""))

  pois.mat <- as.matrix(mod.mcmc)
  pois.dat <- as.data.frame(pois.mat)
  lambda <- pois.dat$lambda
  pois.lambda[faulti,] <- quantile(lambda,c(0.025,0.5,0.975))


  # Gamma renewal process results
  load(paste("../Results/GammaRes/GammaFault-mcmc",faulti,".image",sep=""))

  gam.mat <- as.matrix(mod.mcmc)
  gam.dat <- as.data.frame(gam.mat)
  alpha <- gam.dat$alpha
  beta <- gam.dat$beta
  gam.alpha[faulti,] <- quantile(alpha,c(0.025,0.5,0.975))
  gam.beta[faulti,] <- quantile(beta,c(0.025,0.5,0.975))

 
 # Weibull renewal process results
  load(paste("../Results/WeibRes/WeibFault-mcmc",faulti,".image",sep=""))

  weib.mat <- as.matrix(mod.mcmc)
  weib.dat <- as.data.frame(weib.mat)
  alpha <- weib.dat$alpha
  beta <- weib.dat$beta
  weib.alpha[faulti,] <- quantile(alpha,c(0.025,0.5,0.975))
  weib.beta[faulti,] <- quantile(beta,c(0.025,0.5,0.975))
  

 # BPT renewal process results
  load(paste("../Results/bptRes/bptFault-mcmc",faulti,".image",sep=""))
  load(paste("../Results/bptRes/HPD-toccBPT",faulti,".image",sep=""))
  
  bpt.mat <- as.matrix(mod.mcmc)
  bpt.dat <- as.data.frame(bpt.mat)
  alpha <- bpt.dat$alpha
  mu <- bpt.dat$mu
  bpt.alpha[faulti,] <- quantile(alpha,c(0.025,0.5,0.975))
  bpt.mu[faulti,] <- quantile(mu,c(0.025,0.5,0.975))

  
  # lognormal renewal process results
  load(paste("../Results/lnormRes/lnormFault-mcmc",faulti,".image",sep=""))

  lnorm.mat <- as.matrix(mod.mcmc)
  lnorm.dat <- as.data.frame(lnorm.mat)
  mu <- lnorm.dat$mu
  sig <- lnorm.dat$sig
  lnorm.mu[faulti,] <- quantile(mu,c(0.025,0.5,0.975))
  lnorm.sig[faulti,] <- quantile(sig,c(0.025,0.5,0.975))
 } 

params <- list("pois.lambda"=pois.lambda, "gam.alpha"=gam.alpha, "gam.beta"=gam.beta, 
               "weib.alpha"=weib.alpha, "weib.beta"=weib.beta,
               "lnorm.mu"=lnorm.mu, "lnorm.sig"=lnorm.sig, 
			   "bpt.alpha"=bpt.alpha,"bpt.mu"=bpt.mu)  
  
save(params,file=paste("ParamEstNew.image",sep=""))
 
  









 
  
  