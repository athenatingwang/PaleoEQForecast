## Set the Working Directory to the "/Rcode/RetroForecast" folder
#setwd("/PaleoEQ/Rcode/RetroForecast")


library(rjags)
library(coda)
library(actuar)

a <- dir("../../DataFinal/chronologies_all_final")
nfault <- length(a)

K = 100

Noccur <- NULL
waic <- weights <- matrix(NA,nrow=nfault,ncol=5)


modAveFore <- Fore.Pois <- Fore.gam <- matrix(NA,nrow=nfault,ncol=1500000)
Fore.weib <- Fore.bpt <- Fore.lnorm <- matrix(NA,nrow=nfault,ncol=1500000)

MAbias <- MAvar <- MAmse <- NULL
poisbias <- poisvar <- poismse <- NULL
gambias <- gamvar <- gammse <- NULL
weibbias <- weibvar <- weibmse <- NULL
bptbias <- bptvar <- bptmse <- NULL
lnormbias <- lnormvar <- lnormmse <- NULL


for (faulti in (1:nfault)){
print(faulti)

# Load data
data = read.csv(paste('../../DataFinal/chronologies_all_final/',a[faulti],sep=''), header=FALSE)


# Chronologies have occurrence times in calendar years. 
# m: number of MC samples of the chronologies
# N: Number of earthquakes at the fault
t.occ = data.matrix(data)
m = nrow(t.occ)
N1 <- ncol(t.occ)

# remove the last occurrence time and do NOT use censored time, 
# use the last occurrence time directly to forecast future values
# That is, we don't have a cutoff time
N <- N1-1
N2 <- N-1
t.occ.cens = t.occ[,1:N]


Noccur[faulti] <- N  
inter.cens = t(diff(t(t.occ.cens))) # Transpose, take difference and transpose back


toccurN.sd <- sd(as.vector(t.occ[1:K,N1]))
toccurN.95CI <- as.numeric( diff( quantile(as.vector(t.occ[1:K,N1]),c(0.025,0.975)) )/2 )
if (toccurN.sd <= 1) toccurN.sd <- 1
if (toccurN.95CI <= 1) toccurN.95CI <- 1

## the mean and median of the K MC samples of the last occurrence time that was removed.
toccurN.mean <- mean( as.vector(t.occ[1:K,N1]) )
toccurN.med <- median( as.vector(t.occ[1:K,N1]) )



  # Poisson process results
  load(paste("../../Results/PoisSimRes2/PoisFault-mcmc",faulti,".image",sep=""))

  pois.mat <- as.matrix(mod.mcmc)
  pois.dat <- as.data.frame(pois.mat)
  lambda <- pois.dat$lambda
  n.mc <- length(lambda)
  Fore.Pois[faulti,] <- as.vector(pois.mat[,4:(3+K)])
  Fore.p <- pois.mat[,4:(3+K)]
  
  
  pPois <- plogPois <- array(NA,c(K,N2,length(lambda)))
  ppdPois <- array(NA,c(K,N2))
  pwaicPois <- array(NA,c(K,N2))

  
  for (j in 1:K){ # jth MC sample of the data, K samples in total
    Zj <- eval(parse(text=paste('pois.dat$"Z[',j,']"',sep="")))
    for (i in 1:(N2-1)){ # ith observation of the jth MC sample	
      pPois[j,i,] <- dexp(inter.cens[j,i],lambda*Zj)
      ppdPois[j,i] <- mean(pPois[j,i,])
      plogPois[j,i,] <- dexp(inter.cens[j,i],lambda*Zj,log=TRUE)
      pwaicPois[j,i] <- var(plogPois[j,i,])
    }
    pPois[j,N2,] <- pexp(inter.cens[j,N2],lambda*Zj, lower.tail = FALSE)
    ppdPois[j,N2] <- mean(pPois[j,N2,])
    plogPois[j,N2,] <- pexp(inter.cens[j,N2],lambda*Zj,lower.tail=FALSE,log.p=TRUE)
    pwaicPois[j,N2] <- var(plogPois[j,N2,])
  }
  

  sum(log(ppdPois))
  sum(pwaicPois)

  waicPois <- -2*sum(log(ppdPois)) + 2* sum(pwaicPois)



  # Gamma renewal process results
  load(paste("../../Results/GammaSimRes2/GammaFault-mcmc",faulti,".image",sep=""))

  gam.mat <- as.matrix(mod.mcmc)
  gam.dat <- as.data.frame(gam.mat)
  alpha <- gam.dat$alpha
  beta <- gam.dat$beta
  n.mc <- length(alpha)
  Fore.gam[faulti,] <- as.vector(gam.mat[,6:(5+K)])
  Fore.g <- gam.mat[,6:(5+K)]

  
  pgam <- ploggam <- array(NA,c(K,N2,length(alpha)))
  ppdgam <- array(NA,c(K,N2))
  pwaicgam <- array(NA,c(K,N2))
 
 
  for (j in 1:K){ # jth MC sample of the data, K samples in total
    Yj <- eval(parse(text=paste('gam.dat$"Y[',j,']"',sep="")))
    Zj <- eval(parse(text=paste('gam.dat$"Z[',j,']"',sep="")))
    for (i in 1:(N2-1)){ # ith observation of the jth MC sample	
      pgam[j,i,] <- dgamma(inter.cens[j,i],shape=Zj*alpha,rate=Yj*beta)
      ppdgam[j,i] <- mean(pgam[j,i,])
      ploggam[j,i,] <- dgamma(inter.cens[j,i],shape=Zj*alpha,rate=Yj*beta,log=TRUE)
      pwaicgam[j,i] <- var(ploggam[j,i,])
    }
    pgam[j,N2,] <- pgamma(inter.cens[j,N2],shape=Zj*alpha,rate=Yj*beta,lower.tail=FALSE)
    ppdgam[j,N2] <- mean(pgam[j,N2,])
    ploggam[j,N2,] <- pgamma(inter.cens[j,N2],shape=Zj*alpha,rate=Yj*beta,lower.tail=FALSE,log.p=TRUE)
    pwaicgam[j,N2] <- var(ploggam[j,N2,])

  }
  
  
  sum(log(ppdgam))
  sum(pwaicgam)
  
  waicgam <- -2*sum(log(ppdgam)) + 2* sum(pwaicgam)
  
  
  
  # Weibull renewal process results
  load(paste("../../Results/WeibSimRes2/WeibFault-mcmc",faulti,".image",sep=""))

  weib.mat <- as.matrix(mod.mcmc)
  weib.dat <- as.data.frame(weib.mat)
  alpha <- weib.dat$alpha
  beta <- weib.dat$beta
  n.mc <- length(alpha)
  Fore.weib[faulti,] <- as.vector(weib.mat[,6:(5+K)])
  Fore.w <- weib.mat[,6:(5+K)]
 
  
  pweib <- plogweib <- array(NA,c(K,N2,length(alpha)))
  ppdweib <- array(NA,c(K,N2))
  pwaicweib <- array(NA,c(K,N2))
  

  for (j in 1:K){ # jth MC sample of the data, K samples in total
    Yj <- eval(parse(text=paste('weib.dat$"Y[',j,']"',sep="")))
    Zj <- eval(parse(text=paste('weib.dat$"Z[',j,']"',sep="")))
    shapej <- Zj*alpha
    scalej <- (Yj*beta)^(-1)
    for (i in 1:(N2-1)){ # ith observation of the jth MC sample
      pweib[j,i,] <- dweibull(inter.cens[j,i],shape=shapej,scale=scalej)
      ppdweib[j,i] <- mean(pweib[j,i,])
      plogweib[j,i,] <- dweibull(inter.cens[j,i],shape=shapej,scale=scalej,log=TRUE)
      pwaicweib[j,i] <- var(plogweib[j,i,])
    }
    pweib[j,N2,] <- pweibull(inter.cens[j,N2],shape=shapej,scale=scalej,lower.tail=FALSE)
    ppdweib[j,N2] <- mean(pweib[j,N2,])
    plogweib[j,N2,] <- pweibull(inter.cens[j,N2],shape=shapej,scale=scalej,
                               lower.tail=FALSE,log.p=TRUE)
    pwaicweib[j,N2] <- var(plogweib[j,N2,])    
  }
  
  
  sum(log(ppdweib))
  sum(pwaicweib)
  
  waicweib <- -2*sum(log(ppdweib)) + 2* sum(pwaicweib)
  
  

  # BPT renewal process results
  load(paste("../../Results/bptSimRes2/bptFault-mcmc",faulti,".image",sep=""))
  load(paste("../../Results/bptSimRes2/HPD-toccBPT",faulti,".image",sep=""))
  library(actuar)
  
  bpt.mat <- as.matrix(mod.mcmc)
  bpt.dat <- as.data.frame(bpt.mat)
  alpha <- bpt.dat$alpha
  mu <- bpt.dat$mu
  n.mc <- length(alpha)
  Fore.bpt[faulti,] <- as.vector(bpt.res$t.occ.N)
  Fore.b <- bpt.res$t.occ.N


  pbpt <- plogbpt <- array(NA,c(K,N2,length(alpha)))
  ppdbpt <- array(NA,c(K,N2))
  pwaicbpt <- array(NA,c(K,N2))


  for (j in 1:K){ # jth MC sample of the data, K samples in total
    Yj <- eval(parse(text=paste('bpt.dat$"Y[',j,']"',sep="")))
    Zj <- eval(parse(text=paste('bpt.dat$"Z[',j,']"',sep="")))
    alphaj = alpha * Yj
    muj = mu * Zj
    shapej <- muj/alphaj^2
    
    for (i in 1:(N2-1)){ # ith observation of the jth MC sample
      pbpt[j,i,] <- dinvgauss(inter.cens[j,i], mean=muj, shape = shapej)
      ppdbpt[j,i] <- mean(pbpt[j,i,])
      plogbpt[j,i,] <- dinvgauss(inter.cens[j,i], mean=muj, shape = shapej,log=TRUE)
      pwaicbpt[j,i] <- var(plogbpt[j,i,])
    }
    pbpt[j,N2,] <- pinvgauss(inter.cens[j,N2], mean=muj, shape = shapej,lower.tail=FALSE)
    ppdbpt[j,N2] <- mean(pbpt[j,N2,])
    plogbpt[j,N2,] <- pinvgauss(inter.cens[j,N2], mean=muj, shape = shapej,
                               lower.tail=FALSE,log.p=TRUE)
    pwaicbpt[j,N2] <- var(plogbpt[j,N2,])
  }
  
  
  sum(log(ppdbpt))
  sum(pwaicbpt)
  
  waicbpt <- -2*sum(log(ppdbpt)) + 2* sum(pwaicbpt)
  
  
  
  # lognormal renewal process results
  load(paste("../../Results/lnormSimRes2/lnormFault-mcmc",faulti,".image",sep=""))

  lnorm.mat <- as.matrix(mod.mcmc)
  lnorm.dat <- as.data.frame(lnorm.mat)
  mu <- lnorm.dat$mu
  sig <- lnorm.dat$sig
  n.mc <- length(mu)
  Fore.lnorm[faulti,] <- as.vector(lnorm.mat[,6:(5+K)])
  Fore.l <- lnorm.mat[,6:(5+K)]


  plnorm <- ploglnorm <- array(NA,c(K,N2,length(alpha)))
  ppdlnorm <- array(NA,c(K,N2))
  pwaiclnorm <- array(NA,c(K,N2))
  

  for (j in 1:K){ # jth MC sample of the data, K samples in total
    Yj <- eval(parse(text=paste('lnorm.dat$"Y[',j,']"',sep="")))
    Zj <- eval(parse(text=paste('lnorm.dat$"Z[',j,']"',sep="")))
    meanlj <- Zj+mu
    sdlj <- Yj*sig
    for (i in 1:(N2-1)){ # ith observation of the jth MC sample	
      plnorm[j,i,] <- dlnorm(inter.cens[j,i],meanlj,sdlj)
      ppdlnorm[j,i] <- mean(plnorm[j,i,])
      ploglnorm[j,i,] <- dlnorm(inter.cens[j,i],meanlj,sdlj,log=TRUE)
      pwaiclnorm[j,i] <- var(ploglnorm[j,i,])
    }
    plnorm[j,N2,] <- plnorm(inter.cens[j,N2],meanlj,sdlj,lower.tail=FALSE)
    ppdlnorm[j,N2] <- mean(plnorm[j,N2,])
    ploglnorm[j,N2,] <- plnorm(inter.cens[j,N2],meanlj,sdlj,
                               lower.tail=FALSE,log.p=TRUE)
    pwaiclnorm[j,N2] <- var(ploglnorm[j,N2,])

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
#  modAveFore[faulti,] <- (ind==1)*Fore.Pois[faulti,] + (ind==2)*Fore.gam[faulti,] +
#    (ind==3)*Fore.weib[faulti,] + (ind==4)*Fore.bpt[faulti,] +
#    (ind==5)*Fore.lnorm[faulti,]

  MA.sc <- (ind==1)*Fore.p + (ind==2)*Fore.g + (ind==3)*Fore.w + 
           (ind==4)*Fore.b + (ind==5)*Fore.l

  modAveFore[faulti,] <- as.vector(MA.sc)



poisbias[faulti] <- mean(Fore.Pois[faulti,]) - toccurN.mean
gambias[faulti] <- mean(Fore.gam[faulti,]) - toccurN.mean
weibbias[faulti] <- mean(Fore.weib[faulti,]) - toccurN.mean
bptbias[faulti] <- mean(Fore.bpt[faulti,]) - toccurN.mean
lnormbias[faulti] <- mean(Fore.lnorm[faulti,]) - toccurN.mean
MAbias[faulti] <- mean(modAveFore[faulti,]) - toccurN.mean

poisvar[faulti] <- var(Fore.Pois[faulti,]) 
gamvar[faulti] <- var(Fore.gam[faulti,]) 
weibvar[faulti] <- var(Fore.weib[faulti,]) 
bptvar[faulti] <- var(Fore.bpt[faulti,]) 
lnormvar[faulti] <- var(Fore.lnorm[faulti,])
MAvar[faulti] <- var(modAveFore[faulti,]) 

poismse[faulti] <- poisbias[faulti]^2 + poisvar[faulti]
gammse[faulti] <- gambias[faulti]^2 + gamvar[faulti]
weibmse[faulti] <- weibbias[faulti]^2 + weibvar[faulti]
bptmse[faulti] <- bptbias[faulti]^2 + bptvar[faulti]
lnormmse[faulti] <- lnormbias[faulti]^2 + lnormvar[faulti]
MAmse[faulti] <- MAbias[faulti]^2 + MAvar[faulti]

}

write.csv(cbind(MAbias, MAvar, MAmse, poisbias, poisvar, poismse,
               gambias, gamvar, gammse, weibbias, weibvar, weibmse,
               bptbias, bptvar, bptmse, lnormbias, lnormvar, lnormmse),
               "../../Results/MAforeBiasVarMSE.csv")


ave.foreres <- list("modAveFore"=modAveFore,"Fore.Pois"=Fore.Pois, "Fore.gam"=Fore.gam,
                    "Fore.weib"=Fore.weib,"Fore.bpt"=Fore.bpt,"Fore.lnorm"=Fore.lnorm,
                    "weights"=weights,"waic"=waic) 


save(ave.foreres,file=paste("../../Results/RetroModAveForeRes.image",sep=""))

write.csv(ave.foreres$weights,"../../Results/RetroWAICweights.csv")





gam.alpha <- gam.beta <- weib.alpha <- weib.beta <- matrix(NA,nfault,3)
lnorm.mu <- lnorm.sig <- bpt.alpha <- bpt.mu <- matrix(NA,nfault,3)

for (faulti in 1:nfault){
  # Gamma renewal process results
  load(paste("../../Results/GammaSimRes2/GammaFault-mcmc",faulti,".image",sep=""))

  gam.mat <- as.matrix(mod.mcmc)
  gam.dat <- as.data.frame(gam.mat)
  alpha <- gam.dat$alpha
  beta <- gam.dat$beta
  gam.alpha[faulti,] <- quantile(alpha,c(0.025,0.5,0.975))
  gam.beta[faulti,] <- quantile(beta,c(0.025,0.5,0.975))

 
 # Weibull renewal process results
  load(paste("../../Results/WeibSimRes2/WeibFault-mcmc",faulti,".image",sep=""))

  weib.mat <- as.matrix(mod.mcmc)
  weib.dat <- as.data.frame(weib.mat)
  alpha <- weib.dat$alpha
  beta <- weib.dat$beta
  weib.alpha[faulti,] <- quantile(alpha,c(0.025,0.5,0.975))
  weib.beta[faulti,] <- quantile(beta,c(0.025,0.5,0.975))
  

 # BPT renewal process results
  load(paste("../../Results/bptSimRes2/bptFault-mcmc",faulti,".image",sep=""))
  load(paste("../../Results/bptSimRes2/HPD-toccBPT",faulti,".image",sep=""))
  
  bpt.mat <- as.matrix(mod.mcmc)
  bpt.dat <- as.data.frame(bpt.mat)
  alpha <- bpt.dat$alpha
  mu <- bpt.dat$mu
  bpt.alpha[faulti,] <- quantile(alpha,c(0.025,0.5,0.975))
  bpt.mu[faulti,] <- quantile(mu,c(0.025,0.5,0.975))

  
  # lognormal renewal process results
  load(paste("../../Results/lnormSimRes2/lnormFault-mcmc",faulti,".image",sep=""))

  lnorm.mat <- as.matrix(mod.mcmc)
  lnorm.dat <- as.data.frame(lnorm.mat)
  mu <- lnorm.dat$mu
  sig <- lnorm.dat$sig
  lnorm.mu[faulti,] <- quantile(mu,c(0.025,0.5,0.975))
  lnorm.sig[faulti,] <- quantile(sig,c(0.025,0.5,0.975))
 } 

params <- list("gam.alpha"=gam.alpha, "gam.beta"=gam.beta, 
               "weib.alpha"=weib.alpha, "weib.beta"=weib.beta,
               "lnorm.mu"=lnorm.mu, "lnorm.sig"=lnorm.sig, 
			   "bpt.alpha"=bpt.alpha,"bpt.mu"=bpt.mu)  
  
save(params,file=paste("../../Results/RetroParamEst.image",sep=""))
 





  