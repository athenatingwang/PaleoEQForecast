
###################################################################
###########   Code for all plots and regression results  ##########
### In paper, Bayesian results are used as all posterior samples were 
### used to fit the model in Bayesian method  
###################################################################

library(rjags)
#setwd("B:/Research/PaleoEQ/DataFinal")
setwd("../DataFinal")
a <- dir("chronologies_all_final")
nfault <- length(a)

K = 100


faultname <- c("Alaska PWS Copper", "Alpine Hokuri Ck South Westland", "Awatere East", "Bree",                    
               "Cadell", "Cascadia", "Cascadia Nth", "Cascadia Sth", "Chile Margin",             
               "Cloudy Fault", "Daqingshan Piedmont Hohhot", "Dead Sea Beteiha", "Dead Sea Jordan", "Dead Sea Qatar",            
               "Dead Sea Taybeh", "Dead Sea Yammouneh", "Dunstan", "East Kunlun Kusaihu", "East Kunlun Xidatan",
               "Elashan", "Elsinore Temecula", "Fatigue Wash", "Futugawa", "Garlock El Paso Peaks", 
               "Garlock Twin Lakes", "Gulang Tianqiaogou", "Haiyuan Middle",
               "Hayward Tysons", "Helanshan Eastern Piedmont", "Hope", "Hope Conway", "Hyden", "Irpinia",                 
               "Javon Canyon", "Kiri", "Lachlan", "Lake Edgar", "Langshan Piedmont Xibulong East", "Lenglongling",           
               "Mangatete", "Nankai Trough", "New Guinea", "North Anatolian Cukurcimen",
               "North Anatolian Elmacik", "North Anatolian Gunalan", "North Anatolian Kavakkoy",  
               "North Anatolian Lake Ladik", "North Anatolian Yaylabeli", "Okaya", "Paeroa",                  
               "Pasuruan", "Pihama", "Porters Pass East", "Qilianshan Laohushan", "Rangipo", "Reelfoot", "Rocky Valley",             
               "Rotoitipakau", "San Andreas Big Bend", "San Andreas Burro", "San Andreas Carizzo",       
               "San Andreas Coachella", "San Andreas Mendocino", "San Andreas Mission Ck",     
               "San Andreas Pallet Ck", "San Andreas Pittman", "San Andreas Thousand Palms", 
               "San Andreas Vedanta", "San Andreas Wrightwood", "San Jacinto HogLake",       
               "San Jacinto Mystic Lake", "Serteng Piedmont Wujia", "Snowden", "Solitario Canyon", "Stagecoach Road",          
               "Sumatra Mentawai", "Tanna", "Teton Lakes", "Vernon", "Wairarapa South",          
               "Wairau", "Waitangi", "Wasatch Brigham", "Wasatch Nilphi", "Wasatch Weber",            
               "Wharekuri", "Whirinaki", "Windy Wash", "Wulashan Piedmont Heshunzhuang Botou",
               "Wulashan Piedmont Jinmiaozi Heshuzhuang", "Wutai North Piedmont", "Xorkoli", "Zemuhe") 



Noccur <- NULL
mean.inter <- sd.inter <- NULL

###############################################################
#####  Histogram of inter-arrival times for all faults    #####
###############################################################
#pdf("histInterevent.pdf",paper="special",width=100,height=100)
#par(mfrow=c(10,10))
for (faulti in (1:nfault)){
  data = read.csv(paste('chronologies_all_final/',a[faulti],sep=''), header=FALSE)

# Chronologies have occurrence times in calendar years. 
# m: number of MC samples of the chronologies
# N: Number of earthquakes at the fault
  t.occ = data.matrix(data)
  m = nrow(t.occ)
  N1 <- ncol(t.occ)
  
## Interarrival times mean and variance
  diffoccur <- apply(t.occ[1:100,],1,diff)
  mean.inter <- append(mean.inter,mean(diffoccur))
  sd.inter <- append(sd.inter,sd(diffoccur))
#  hist(diffoccur,cex.axis=4,main=faultname[faulti],cex=4,xlim=c(0,max(diffoccur)))
  
# remove the last occurrence time and do NOT use censored time, 
# use the last occurrence time directly to forecast future values
# That is, we don't have a cutoff time
  N <- N1
  N2 <- N-1
  Noccur[faulti] <- N  
}
#dev.off()



###############################################################
#####################   Read Results   ########################
###############################################################
# setwd("B:/Research/PaleoEQ/ResultsFinal")
setwd("../ResultsFinal")
load(paste("ModAveProbResNew.image",sep=""))
load(paste("ModAveLogProbResNew.image",sep=""))
load(paste("ModAveForeResNew.image",sep=""))

weights <- read.csv("WAICweightsNew.csv")

## Load the file containing the quantiles of the paramter estimates
load(paste("ParamEstNew.image",sep=""))


## Save the quantiles of each forecast of natual log probability of an earthquake in the next 50 years on each fault.
lp50Pois.quant <- lp50gam.quant <- lp50weib.quant <- lp50bpt.quant <- lp50lnorm.quant <- lp50modAve.quant <- matrix(NA,nfault,3)
for (faulti in 1:nfault){
  lp50Pois.quant[faulti,] <- quantile(as.vector(ave.logprobres$lp50Pois[faulti,,]),c(0.025,0.5,0.975))
  lp50gam.quant[faulti,] <- quantile(as.vector(ave.logprobres$lp50gam[faulti,,]),c(0.025,0.5,0.975))
  lp50weib.quant[faulti,] <- quantile(as.vector(ave.logprobres$lp50weib[faulti,,]),c(0.025,0.5,0.975))
  lp50bpt.quant[faulti,] <- quantile(as.vector(ave.logprobres$lp50bpt[faulti,,]),c(0.025,0.5,0.975))
  lp50lnorm.quant[faulti,] <- quantile(as.vector(ave.logprobres$lp50lnorm[faulti,,]),c(0.025,0.5,0.975))
  lp50modAve.quant[faulti,] <- quantile(as.vector(ave.logprobres$lp50modAve[faulti,,]),c(0.025,0.5,0.975))
}
Forelp50quantile <- list("lp50Pois.quant"=lp50Pois.quant,"lp50gam.quant"=lp50gam.quant,
          "lp50weib.quant"=lp50weib.quant,"lp50bpt.quant"=lp50bpt.quant,
          "lp50lnorm.quant"=lp50lnorm.quant,"lp50modAve.quant"=lp50modAve.quant)

#save(Forelp50quantile ,file="Forecastlogprob50yr.image")



p50Pois.quant <- p50gam.quant <- p50weib.quant <- p50bpt.quant <- p50lnorm.quant <- p50modAve.quant <- matrix(NA,nfault,3)
for (faulti in 1:nfault){
  p50Pois.quant[faulti,] <- quantile(as.vector(ave.probres$p50Pois[faulti,,]),c(0.025,0.5,0.975))
  p50gam.quant[faulti,] <- quantile(as.vector(ave.probres$p50gam[faulti,,]),c(0.025,0.5,0.975))
  p50weib.quant[faulti,] <- quantile(as.vector(ave.probres$p50weib[faulti,,]),c(0.025,0.5,0.975))
  p50bpt.quant[faulti,] <- quantile(as.vector(ave.probres$p50bpt[faulti,,]),c(0.025,0.5,0.975))
  p50lnorm.quant[faulti,] <- quantile(as.vector(ave.probres$p50lnorm[faulti,,]),c(0.025,0.5,0.975))
  p50modAve.quant[faulti,] <- quantile(as.vector(ave.probres$p50modAve[faulti,,]),c(0.025,0.5,0.975))
}
Forep50quantile <- list("p50Pois.quant"=p50Pois.quant,"p50gam.quant"=p50gam.quant,
          "p50weib.quant"=p50weib.quant,"p50bpt.quant"=p50bpt.quant,
          "p50lnorm.quant"=p50lnorm.quant,"p50modAve.quant"=p50modAve.quant)

#save(Forep50quantile ,file="Forecastprob50yr.image")


xxx <- formatC(as.matrix(weights[,2:6]), format = "f", digits = 4)
#colnames(xxx)=c("wtPois","wtgam","wtweib","wtbpt","wtlnorm")
# write.csv(as.data.frame(cbind(faultname,xxx)),"WeightsFaultname.csv")

#######################################################################################


cbind(faultname, formatC(p50modAve.quant, format = "f", digits = 4),mean.inter,sd.inter)
cbind(faultname, formatC(p50weib.quant, format = "f", digits = 4),mean.inter,sd.inter)
cbind(faultname, formatC(p50Pois.quant, format = "f", digits = 4),mean.inter,sd.inter)
cbind(faultname, formatC(p50gam.quant, format = "f", digits = 4),mean.inter,sd.inter)
cbind(faultname, formatC(p50lnorm.quant, format = "f", digits = 4),mean.inter,sd.inter)
cbind(faultname, formatC(p50bpt.quant, format = "f", digits = 4),mean.inter,sd.inter)



###
cut <- 0.95

ind <- 1:nfault
ind.Pois <- ind[weights$V1>=cut]  
ind.gam <- ind[weights$V2>=cut]  
ind.weib <- ind[weights$V3>=cut]  
ind.bpt <- ind[weights$V4>=cut]  
ind.lnorm <- ind[weights$V5>=cut]  
ind.ma <- ind[-sort(c(ind.Pois,ind.gam,ind.weib,ind.bpt,ind.lnorm))]  

weights[c(12,14,22,24,51,56,58,86),]



length(ind.Pois)
length(ind.gam)
length(ind.weib)
length(ind.bpt)
length(ind.lnorm)


faultname[ind.weib[params$weib.alpha[ind.weib,3]<=1]]
faultname[ind.gam[params$gam.alpha[ind.gam,3]<=1]]

faultname[ind.ma[params$gam.alpha[ind.ma,3]<=1]]
faultname[ind.ma[params$weib.alpha[ind.ma,3]<=1]]

faultname[ind.Pois]


######################################################################
#######  Figure 1 bar plots
######################################################################
x1 <- matrix(0,nrow=5,ncol=31)
temp <- as.data.frame(table(Noccur[ind.ma]))
x1[1,as.numeric(as.character(temp[,1]))-4] <- temp[,2]
temp <- as.data.frame(table(Noccur[ind.gam]))
x1[2,as.numeric(as.character(temp[,1]))-4] <- temp[,2]
temp <- as.data.frame(table(Noccur[ind.weib]))
x1[3,as.numeric(as.character(temp[,1]))-4] <- temp[,2]
temp <- as.data.frame(table(Noccur[ind.bpt]))
x1[4,as.numeric(as.character(temp[,1]))-4] <- temp[,2]
temp <- as.data.frame(table(Noccur[ind.lnorm]))
x1[5,as.numeric(as.character(temp[,1]))-4] <- temp[,2]
modelapp <- c("MA","Gamma","Weibull","BPT","Lognormal") 
colors = c("gray","purple","green","blue","brown")

x2 <- x1/matrix(rep(colSums(x1),5),nrow=5,byrow=T)


add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}
mods <- c("Poisson","Gamma","Weibull","BPT","Lognormal") 
colors2 = c("red","purple","green","blue","brown")

allmodels <- c("MA","Poisson","Gamma","Weibull","BPT","Lognormal") 
colors.lgd = c("gray","red","purple","green","blue","brown")

postscript("BarplotPreferredModel-MdWt.eps",paper="special",width=5*1.6*2*0.8,height=5*0.9*2,onefile = TRUE, horizontal = FALSE)
layout( matrix(c(1,1,2,3), 2,2,byrow=TRUE) )
par(mar=c(4.5,4.5,4,1))
barplot(as.matrix(t(weights[,2:6])),main="(a)",names.arg = seq(1,nfault,1),xlab="Fault No",
        ylab="Weights",col=colors2,cex=1.5,cex.lab=1.5,cex.axis=1.5)
box()
#add_legend("top", legend=mods, cex = 1.5, fill = colors,horiz=TRUE, bty='n')

par(mar=c(4.5,4.5,4,1))
barplot(x1,main="(b)",names.arg = seq(5,35,1),xlab="Number of events",
        ylab="Frequency",col=colors,cex=1.5,cex.lab=1.5,cex.axis=1.5)
box()
legend("topright", legend=allmodels, cex = 1.5, fill = colors.lgd)

par(mar=c(4.5,4.5,4,1))
barplot(x2,main="(c)",names.arg = seq(5,35,1),xlab="Number of events",
        ylab="Proportion",col=colors,cex=1.5,cex.lab=1.5,cex.axis=1.5)
#legend("topright", legend=modelapp, cex = 1.5, fill = colors)
box()
dev.off()
######################################################################



ind.periodic <- ind[ params$weib.alpha[,1]>=1 &
                     params$gam.alpha[,1]>=1 ]
 
ind.cluster <- ind[ params$weib.alpha[,3] < 1 &
                     params$gam.alpha[,3] < 1 ]
 
## Near Poisson behaviour
faultname[ind[-c(ind.periodic,ind.cluster)]]
  
params$gam.alpha[ind.cluster,]  
params$weib.alpha[ind.cluster,]
params$bpt.alpha[ind.cluster,]



###########################################################################
## Get the single best model, and forecast quantiles from single best model
###########################################################################
maxweight <- NULL
for (i in 1:nfault){
  maxweight[i] <- which.max(weights[i,2:6])
}


bmforequant <- Forep50quantile$p50Pois.quant
for (i in 1:nfault){
  if (maxweight[i]==1) bmforequant[i,] <- Forep50quantile$p50Pois.quant[i,]
  if (maxweight[i]==2) bmforequant[i,] <- Forep50quantile$p50gam.quant[i,]
  if (maxweight[i]==3) bmforequant[i,] <- Forep50quantile$p50weib.quant[i,]
  if (maxweight[i]==4) bmforequant[i,] <- Forep50quantile$p50bpt.quant[i,]
  if (maxweight[i]==5) bmforequant[i,] <- Forep50quantile$p50lnorm.quant[i,]
}

bestmod <- NULL
for (i in 1:nfault){
  if (maxweight[i]==1) bestmod[i] <- "P"
  if (maxweight[i]==2) bestmod[i] <- "G"
  if (maxweight[i]==3) bestmod[i] <- "W"
  if (maxweight[i]==4) bestmod[i] <- "B"
  if (maxweight[i]==5) bestmod[i] <- "L"
}

########################################################################
############  Parameter estimates for all San Andreas segments  ########
########################################################################
params$gam.alpha[59:69,]
params$weib.alpha[59:69,]
params$bpt.alpha[59:69,]
Noccur[59:69]
bestmod[59:69]




#########################################################################
###   The faults that have the Weibull renewal process as the single-
###   best model have on average about 1.45 (90% CI (1.04,2.08)) times 
###   more large earthquakes than those for which the BPT renewal process 
###   is the single-best model
#########################################################################

mod.best <- rep(NA,nfault)
mod.best[ind.Pois] <- "P"
mod.best[ind.gam] <- "G"
mod.best[ind.weib] <- "W"
mod.best[ind.bpt] <- "B"
mod.best[ind.lnorm] <- "L"
mod.best[-sort(c(ind.Pois,ind.gam,ind.weib,ind.bpt,ind.lnorm))] <- "MA"
mod.best.no <- as.numeric(as.factor(mod.best))

write.csv(rbind(seq(5,35,1),BPT,Gam,Lnorm,MA,Pois,Weib),"NoccurBestMod.csv")

########################################################################
#########  Relation between number of events and best model
########################################################################

##### Bayesian method
tmp <- "
model{
  for (j in 1:nfault){ # jth fault
    Noccur[j] ~ dnegbin(p[j],size)
    p[j] <- size/(size+mu[j])
    log(mu[j]) <- a + b[X[j]]
  }
  a ~ dnorm(0,1.0E-06)
  b[1] <- 0
  for (i in 2:5){
    b[i] ~ dnorm(0,1.0E-06)
  }
  size ~ dunif(0.001,1000)
  theta <- pow(1/mean(p),2)
  scaleparam <- mean((1-p)/p) 
}
"

cat(tmp, file = "NoccurBestMod.Negbin.jags")

library(R2jags)
X <- mod.best.no
jagsdata <- list("Noccur", "X", "nfault")
params <- c("a","b[2]","b[3]","b[4]","b[5]","size","theta","scaleparam")

system.time(
  NB.Noccur <- jags(data=jagsdata, inits=NULL,
                    parameters.to.save = params,n.chain=3,
                    n.iter=20000, n.thin=10, n.burnin=10000,
                    model.file = 'NoccurBestMod.Negbin.jags')
)

# Convert to an MCMC object
mod.mcmc <- as.mcmc(NB.Noccur)
plot(mod.mcmc,trace=TRUE)
gelman.diag(mod.mcmc)
summary(mod.mcmc)

mod.mat <- as.matrix(mod.mcmc)
mod.dat <- as.data.frame(mod.mat)

## Quantile of coefficient exp(b)
expb.quant <- matrix(NA,nrow=4,ncol=5)
for (m in 2:5){
  bm <- eval(parse(text=paste('mod.dat$"b[',m,']"',sep="")))
  expb.quant[m-1,] <- quantile(exp(bm),c(0.025,0.05,0.5,0.95,0.975))
}
expb.quant

                   [,1]      [,2]     [,3]     [,4]     [,5]
b2 Gamm  [1,] 0.8439750 0.8979994 1.235469 1.758345 1.852133
b3 Lnorm [2,] 0.5694068 0.6254751 1.099618 1.895289 2.092137
b4 MA    [3,] 0.7101419 0.7608302 1.053136 1.470018 1.568064
b5 Weib  [4,] 1.0336945 1.0878492 1.449139 1.945754 2.064496 *******



#####  Frequentist method
dat <- as.data.frame(cbind(Noccur,mod.best))
dat$Noccur <- as.numeric(dat$Noccur)
dat$mod.best <- as.factor(dat$mod.best)

a <- glm(Noccur~mod.best,family="quasipoisson")
summary(a)
anova(a,test="Chisq")

est <- exp(coef(a))
cf95 <- confint(a,level=0.95)
cf90 <- confint(a,level=0.9)

cfex95 <- exp(cf95)
cfex90 <- exp(cf90)

cbind(est,cfex95,cfex90)


#postscript("BestModelNoEventsGLM.eps",paper="special",width=5*1.6*2,height=5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,2))
par(mar=c(4.5,4.5,1,1))
rsd <- rstandard(a)
ks.test(rsd,"pnorm")
lp<-predict(a)
plot(lp,rsd,xlab="Predited values",ylab="Residuals")
abline(h=0)
qqnorm(rsd,xlab="...",ylab="...",main="")
abline(0,1)  
#dev.off()
#########################################################################





# Section 2

########################################################################
###############    Table 1
####  Print fault name and MA forecast and Best model forecast
########################################################################
cbind(faultname,Noccur,formatC(as.matrix(Forep50quantile$p50modAve.quant), format = "f", digits = 4),
       bestmod,formatC(as.matrix(bmforequant), format = "f", digits = 4))





########################################################################
## Difference between averaged and best model 50 year forecast probability
########################################################################
which(abs(Forep50quantile$p50modAve.quant[,3]-bmforequant[,3])>=0.001)
#[1] 21 42 56 61 67 81 91
which(abs(Forep50quantile$p50modAve.quant[,1]-bmforequant[,1])>=0.001)
#[1] 12 14 30 61 64 67 81
########################################################################



########################################################################
######  Faults with large probability of eq occurrence and very low 
######  probability of eq occurrence in the next 50 years
########################################################################
faultname[which(Forep50quantile$p50modAve.quant[,2]>=0.5)]
#[1] "San Andreas Carizzo"    "San Andreas Wrightwood" "Sumatra Mentawai"
faultname[which(Forep50quantile$p50modAve.quant[,3]<=0.00001)]
#[1] "East Kunlun Kusaihu"        "Futugawa"                  
#[3] "Helanshan Eastern Piedmont" "Qilianshan Laohushan" 
########################################################################






#########################################################################
######   On average, faults along a plate boundary are about 8.37 (90% CI 
######   (5.37,13.11)) times more likely to have a large earthquake in 
######   the next 50 years than intraplate faults.
#########################################################################
###  Load meta data on tectonic settings and faulting style.
# setwd("B:/Research/PaleoEQ/DataFinal")
meta <- read.csv("metadata_summary.csv")
attach(meta)
names(meta)
Tec_region <- meta$Tectonic_region
Fault_style <- meta$Faulting_style

tecreg <- as.numeric(as.factor(Tec_region))
unique(tecreg)
unique(Tec_region)

#> unique(tecreg)
#[1] 7 5 6 3 1 4 2
#> unique(Tectonic_region)
#[1] "Subduction"             "Plate_boundary_master"  "Plate_boundary_network"
#[4] "Intraplate_noncratonic" "Active_intraplate"      "Near_plate_boundary"
#[7] "Intraplate_cratonic"

Tec_region[Tec_region=="Subduction"] <- "S"
Tec_region[Tec_region=="Plate_boundary_master"] <- "B"
Tec_region[Tec_region=="Plate_boundary_network"] <- "B"
Tec_region[Tec_region=="Intraplate_noncratonic"] <- "I"
Tec_region[Tec_region=="Active_intraplate"] <- "AI"
Tec_region[Tec_region=="Near_plate_boundary"] <- "NB"
Tec_region[Tec_region=="Intraplate_cratonic"] <- "I"

tecreg12fc <- Tec_region
tecreg12fc[Tec_region=="S"] <- "PB"
tecreg12fc[Tec_region=="B"] <- "PB"
tecreg12fc[Tec_region=="B"] <- "PB"
tecreg12fc[Tec_region=="I"] <- "IP"
tecreg12fc[Tec_region=="AI"] <- "IP"
tecreg12fc[Tec_region=="NB"] <- "IP"
tecreg12fc[Tec_region=="I"] <- "IP"


### Load probability of an event in the next 50 years
# setwd("B:/Research/PaleoEQ/ResultsFinal")
# load(paste("ModAveProbResNew.image",sep=""))
nmcmc <- length(as.vector(ave.probres$p50modAve[1,,]))
p50modAve <- matrix(NA,nfault,nmcmc)
for (faulti in 1:nfault){
  p50modAve[faulti,] <- as.vector(ave.probres$p50modAve[faulti,,])+0.000001
}
#ks.test(log(p50modAve[faulti,]),"pnorm",mean=mean(log(p50modAve[3,])),sd=sd(log(p50modAve[faulti,])))



############################################################################
####  NOT USED IN PAPER. SEE LOGNORMAL LATER THAT HAS SMALLER WAIC
############################################################################
### Bayesian Method Gamma
tmp1 <- "
model{
  for (i in 1:nmcmc){
    for (j in 1:nfault){ # jth fault
     p50modAve[j,i] ~ dgamma(shape,shape/exp(logmu[j,i]))
     logmu[j,i] <- a[i] + b[tecreg12[j],i]
   }
   a[i] <- mua + epsa[i]
   epsa[i] ~ dnorm(0,1/asig/asig)
   b[2,i] <- mub + epsb[i]
   epsb[i] ~ dnorm(0,1/bsig/bsig)
   b[1,i] <- 0
   expb[i] <- exp(b[2,i])
  }
  mua ~ dnorm(0,1.0E-06)
  asig ~ dt(0,0.04,3)T(0,)
  bsig ~ dt(0,0.04,3)T(0,)
  mub ~ dnorm(0,1.0E-06)
  shape ~ dunif(0,100)
}
"

cat(tmp1, file = "Pro50tecregGam.jags")

library(R2jags)
nmcmc <- 1000
tecreg12 <- as.numeric(as.factor(tecreg12fc))
unique(tecreg12fc)
unique(tecreg12)
jagsdata <- list("p50modAve", "tecreg12", "nfault","nmcmc")
params <- c("mua","mub","shape","asig","bsig","a","b")


system.time(
  gam.prob <- jags(data=jagsdata, inits=NULL,
                     parameters.to.save = params,n.chain=3,
                     n.iter=25000, n.thin=10, n.burnin=5000,
                     model.file = 'Pro50tecregGam.jags')
)

# Convert to an MCMC object
mod.mcmc <- as.mcmc(gam.prob)

save(mod.mcmc,file=paste("Pro50tecregRegressGam-mcmc.image",sep=""))


#plot(mod.mcmc,trace=TRUE)
#gelman.diag(mod.mcmc)
#summary(mod.mcmc)

mod.mat <- as.matrix(mod.mcmc)
mod.dat <- as.data.frame(mod.mat)

## Quantile of coefficient exp(b)
expball <- c()
for (i in 1:nmcmc){
  bm <- eval(parse(text=paste('mod.dat$"b[2,',i,']"',sep="")))
  expball <- append(expball,exp(bm))
}
expb.quant <- quantile(expball,c(0.025,0.05,0.5,0.95,0.975))
expb.quant


bm <- eval(parse(text=paste('mod.dat$"mub"',sep="")))
expb.quant <- quantile(exp(bm),c(0.025,0.05,0.5,0.95,0.975))
expb.quant


###  WAIC
get(load("Pro50tecregRegressGam-mcmc.image"))
mod.mat <- as.matrix(mod.mcmc)
mod.dat <- as.data.frame(mod.mat)

gamshape <- eval(parse(text=paste('mod.dat$shape',sep="")))
b1 <- rep(0,length(gamshape))
logmu <- pgam <- array(NA,c(nfault,nmcmc,length(gamshape)))
pwaicgam <- array(NA,c(nfault,nmcmc))
for (i in 1:nmcmc){
  b2 <- eval(parse(text=paste('mod.dat$"b[2,',i,']"',sep="")))
  bm <- rbind(b1,b2)
  am <- eval(parse(text=paste('mod.dat$"a[',i,']"',sep="")))
  for (j in 1:nfault){ # jth fault
    logmu[j,i,] <- am + bm[tecreg12[j],i]
    pgam[j,i,] <- dgamma(p50modAve[j,i],gamshape,gamshape/exp(logmu[j,i,]))
    pwaicgam[j,i] <- var(dgamma(p50modAve[j,i],gamshape,gamshape/exp(logmu[j,i,]),log = T))
  }
}
sum(log(apply(pgam,3,mean)))
sum(pwaicgam)

waicgam <- -2*sum(log(apply(pgam,3,mean))) + 2* sum(pwaicgam)
waicgam
#-67231.85
############################################################################
####  NOT USED IN PAPER. SEE LOGNORMAL LATER THAT HAS SMALLER WAIC
############################################################################





############################################################################
####  Paper used lognormal model as its WAIC is much smaller than gamma one
############################################################################
##### Bayesian method lognormal
tmp1 <- "
model{
  for (i in 1:nmcmc){
    for (j in 1:nfault){ # jth fault
     p50modAve[j,i] ~ dlnorm(mu[j,i],tau)
     mu[j,i] <- a[i] + b[tecreg12[j],i]
   }
   a[i] <- mua + epsa[i]
   epsa[i] ~ dnorm(0,1/asig/asig)
   b[2,i] <- mub + epsb[i]
   epsb[i] ~ dnorm(0,1/bsig/bsig)
   b[1,i] <- 0
   expb[i] <- exp(b[2,i])
  }
  mua ~ dnorm(0,1.0E-06)
  asig ~ dt(0,0.04,3)T(0,)
  bsig ~ dt(0,0.04,3)T(0,)
  mub ~ dnorm(0,1.0E-06)
  tau ~ dt(0,0.04,3)T(0,)
}
"

cat(tmp1, file = "Pro50tecreg.jags")

library(R2jags)
nmcmc <- 1000
tecreg12 <- as.numeric(as.factor(tecreg12fc))
unique(tecreg12fc)
unique(tecreg12)
jagsdata <- list("p50modAve", "tecreg12", "nfault","nmcmc")
params <- c("mua","mub","tau","asig","bsig","a","b")


system.time(
  lnorm.prob <- jags(data=jagsdata, inits=NULL,
                     parameters.to.save = params,n.chain=3,
                     n.iter=25000, n.thin=10, n.burnin=5000,
                     model.file = 'Pro50tecreg.jags')
)

# Convert to an MCMC object
mod.mcmc <- as.mcmc(lnorm.prob)

save(mod.mcmc,file=paste("Pro50tecregRegress-mcmc.image",sep=""))


plot(mod.mcmc,trace=TRUE)
gelman.diag(mod.mcmc)
summary(mod.mcmc)

mod.mat <- as.matrix(mod.mcmc)
mod.dat <- as.data.frame(mod.mat)

## Quantile of coefficient exp(b)
expball <- c()
for (i in 1:nmcmc){
  bm <- eval(parse(text=paste('mod.dat$"b[2,',i,']"',sep="")))
  expball <- append(expball,exp(bm))
}
expb.quant <- quantile(expball,c(0.025,0.05,0.5,0.95,0.975))
expb.quant
# 25.44072 25.64546 26.62353 27.64157 27.85569

bm <- eval(parse(text=paste('mod.dat$"mub"',sep="")))
expb.quant <- quantile(exp(bm),c(0.025,0.05,0.5,0.95,0.975))
expb.quant



###  WAIC
get(load("Pro50tecregRegress-mcmc.image"))
mod.mat <- as.matrix(mod.mcmc)
mod.dat <- as.data.frame(mod.mat)

lnormsig <- 1/eval(parse(text=paste('mod.dat$tau',sep="")))
b1 <- rep(0,length(lnormsig))
mu <- plnm <- array(NA,c(nfault,nmcmc,length(lnormsig)))
pwaiclnm <- array(NA,c(nfault,nmcmc))
for (i in 1:nmcmc){
  b2 <- eval(parse(text=paste('mod.dat$"b[2,',i,']"',sep="")))
  bm <- rbind(b1,b2)
  am <- eval(parse(text=paste('mod.dat$"a[',i,']"',sep="")))
  for (j in 1:nfault){ # jth fault
    mu[j,i,] <- am + bm[tecreg12[j],i]
    plnm[j,i,] <- dlnorm(p50modAve[j,i],mu[j,i,],lnormsig)
    pwaiclnm[j,i] <- var(dlnorm(p50modAve[j,i],mu[j,i,],lnormsig,log = T))
  }
}
sum(log(apply(plnm,3,mean)))
sum(pwaiclnm)

waiclnm <- -2*sum(log(apply(plnm,3,mean))) + 2* sum(pwaiclnm)
###  -86587.05
#####################################################################





###################################################################
####################        NOT  USED          ####################
###################################################################
tmp1 <- "
model{
  for (i in 1:nmcmc){
    for (j in 1:nfault){ # jth fault
     p50modAve[j,i] ~ dlnorm(mu[j,i],tau)
     mu[j,i] <- a + b[tecreg12[j]]
   }
  }
  a ~ dnorm(0,1.0E-06)
  b[1] <- 0
  b[2] ~ dnorm(0,1.0E-06)
  tau ~ dt(0,0.04,3)T(0,)
}
"

cat(tmp1, file = "Pro50tecregsimp.jags")

library(R2jags)
nmcmc <- 1000
tecreg12 <- as.numeric(as.factor(tecreg12fc))
unique(tecreg12fc)
unique(tecreg12)
jagsdata <- list("p50modAve", "tecreg12", "nfault","nmcmc")
params <- c("a","b[2]","tau")

system.time(
  lnorm.prob <- jags(data=jagsdata, inits=NULL,
                     parameters.to.save = params,n.chain=3,
                     n.iter=15000, n.thin=10, n.burnin=5000,
                     model.file = 'Pro50tecregsimp.jags')
)


# Convert to an MCMC object
mod.mcmc <- as.mcmc(lnorm.prob)
plot(mod.mcmc,trace=TRUE)
gelman.diag(mod.mcmc)
summary(mod.mcmc)

save(mod.mcmc,file=paste("Pro50tecregRegressSimp-mcmc.image",sep=""))

mod.mat <- as.matrix(mod.mcmc)
mod.dat <- as.data.frame(mod.mat)

## Quantile of coefficient exp(b)
bm <- eval(parse(text=paste('mod.dat$"b[',2,']"',sep="")))
expb.quant <- quantile(exp(bm),c(0.025,0.05,0.5,0.95,0.975))
expb.quant
###################################################################
####################        NOT  USED          ####################
###################################################################




#########     Frequentist method
##  Relation between median probability of an event in 50 years and tectonic region
##  lognormal regression better than gamma regression
p50ma <- Forep50quantile$p50modAve.quant[,2]
p50ma[p50ma<=1e-5] <- 1e-5

dat <- as.data.frame(cbind(p50ma,tecreg12fc))

dat$p50ma<- as.numeric(dat$p50ma)
dat$tecreg12fc <- as.factor(dat$tecreg12fc )
dat$tecreg12fc <- relevel(dat$tecreg12fc,ref="IP")


b <- glm(log(p50ma)~tecreg12fc,data=dat)
summary(b)

est <- exp(coef(b))
cf95 <- confint(b,level=0.95)
cf90 <- confint(b,level=0.9)

cfex95 <- exp(cf95)
cfex90 <- exp(cf90)

cbind(est,cfex95,cfex90)

rsd <- rstandard(b)
ks.test(rsd,"pnorm")
########################################################################






#########################################################################
#############               Figure 2                       ##############
#########################################################################

## Median of model averaged forecast
fmed <- Forelp50quantile$lp50modAve.quant[,2]/2.302585 ## (convert from natual log to log10)
fmed[fmed<=-10000] <- -10000 ## Give -Inf a finite random value
## Problog10 <- fmed
## mapProbdat <- as.data.frame(cbind(Longitude,Latitude,Problog10))
## head(mapProbdat)
## write.csv(mapProbdat,"Prob50yearMetadata.csv")

## 95% credible intervals
f0.025 <- Forelp50quantile$lp50modAve.quant[,1]/2.302585
f0.025[f0.025<=-10000] <- -10000
f0.975 <- Forelp50quantile$lp50modAve.quant[,3]/2.302585
f0.975[f0.975<=-10000] <- -10000


col.med <- rep(NA,nfault)
col.025 <- rep(NA,nfault)
col.975 <- rep(NA,nfault)

library(wesanderson)
ngroup = 8
paltem <- wes_palette("Zissou1", ngroup+4, type = "continuous")
pal <- paltem[c(1:4,7,9,10,12)]

#paltem <- rev(rainbow(n=ngroup*2, start=0, end=4/6))
#pal <- paltem[seq(1,ngroup*2,2)]

mgrid <- cut(fmed,c(-100000,seq(-5,-1,1),log10(0.25),log10(0.5),log10(0.75),0),include.lowest = FALSE, right = TRUE)
col.med <- as.numeric(mgrid)
for (i in 1:ngroup){
  col.med[col.med==i] <- pal[i]
}  


mgrid <- cut(f0.025,c(-100000,seq(-5,-1,1),log10(0.25),log10(0.5),log10(0.75),0),include.lowest = FALSE, right = TRUE)
col.025 <- as.numeric(mgrid)
for (i in 1:ngroup){
  col.025[col.025==i] <- pal[i]
}  


mgrid <- cut(f0.975,c(-100000,seq(-5,-1,1),log10(0.25),log10(0.5),log10(0.75),0),include.lowest = FALSE, right = TRUE)
col.975 <- as.numeric(mgrid)
for (i in 1:ngroup){
  col.975[col.975==i] <- pal[i]
}  

Prob <- c(expression(10^-5),expression(10^-4),expression(10^-3),
        expression(10^-2),expression(10^-1),0.25,0.5,0.75) #,1


library(maps)


postscript("Faultmapp50yr.eps",paper="special",width=5*1.6,height=5*2*0.9,onefile = TRUE, horizontal = FALSE)
layout( matrix(c(1,1,2,3), 2,2,byrow=TRUE),  widths=c(5,4), 
  heights=c(2,2) )
par(mar=c(4.5,4.5,1,1))
map('world', resolution=1, fill = T, 
    col = "gray", border="darkgrey")
map.axes(cex.axis=1.5)
points(Longitude,Latitude,pch=16,col=col.med,cex=1)
legend("bottomleft",legend=Prob,col=pal,pch=16,cex=1.5,bty="n")
text(-170,75,"(a)",cex=1.5)

par(mar=c(4.5,4.5,1,1))
map("world", resolution=1, fill = T, 
    col = "gray", border="darkgrey", xlim=c(-123.5,-115),ylim=c(32.5,39))
map.axes(cex.axis=1.5)
points(Longitude,Latitude,pch=16,col=col.med,cex=1.7)
#legend("bottomleft",legend=Prob,col=pal,pch=16,cex=1.5,bty="n")
text(-123,38.5,"(b)",cex=1.5)

par(mar=c(4.5,4.5,1,1))
map("world", resolution=1, fill = T, 
    col = "gray", border="darkgrey",xlim=c(160,180),ylim=c(-49,-34))
map.axes(cex.axis=1.5)
points(Longitude,Latitude,pch=16,col=col.med,cex=1.7)
#legend("bottomleft",legend=Prob,col=pal,pch=16,cex=1.5,bty="n")
text(162,-35,"(c)",cex=1.5)
dev.off()
##################################################################






#########################################################################
#############               Figure 3                       ##############
#########################################################################

##############################################
### posterior densities of MA forecasts for all faults 
##############################################
fore.range <- t(apply(ave.foreres$modAveFore,1,range))
fore.quant <- t(apply(ave.foreres$modAveFore,1,quantile,probs =0.95))

ind <- 1:nfault
ind1 <- ind[fore.quant<=3000]
ind2 <- ind[fore.quant>3000 & fore.quant<=5000]
ind3 <- ind[fore.quant>5000]

n1 <- length(ind1)
n2 <- length(ind2)
n3 <- length(ind3)

faultname <- c("Alaska", "Alpine", "Awatere", "Bree","Cadell", "Cascadia", "CascaNth", "CascaSth", "Chile",             
               "Cloudy", "Daqing", "DS Bet", "DS Jor", "DS Qat","DS Tay", "DS Yam", "Dunstan", "Kusaihu", "Xidatan",
               "Elashan", "Elsinore", "Fatigue", "Futugawa", "Garlock E","Garlock T", "Gulang", "Haiyuan",
               "Hayward", "Helanshan", "Hope", "Hope Con", "Hyden", "Irpinia",                 
               "Javon", "Kiri", "Lachlan", "LakeEdgar", "Langshan", "LengLL",           
               "Mangatete", "Nankai", "NewGuin", "NA Cuk",
               "NA Elm", "NA Gun", "NA Kav", "NA Lak", "NA Yay", "Okaya", "Paeroa",                  
               "Pasuruan", "Pihama", "Porters", "Qilian", "Rangipo", "Reelfoot", "Rocky",             
               "Rotoit", "SA Big", "SA Bur", "SA Car","SA Coa", "SA Men", "SA Mis",     
               "SA Pal", "SA Pit", "SA Tho","SA Ved", "SA Wri", "SJ Hog",       
               "SJ Mys", "Serteng", "Snowden", "Solitario", "Stagecoach",          
               "Sumatra", "Tanna", "Teton", "Vernon", "Wairarapa",          
               "Wairau", "Waitangi", "WasatchB", "WasatchN", "WasatchW",            
               "Wharekuri", "Whirinaki", "Windy", "Wula-Botou",
               "Wula-Jinm", "Wutai", "Xorkoli", "Zemuhe")


library(coda)
library(R2jags)
library(lattice)

library("bayesplot")
library("ggplot2")
library("rstanarm")
library(grid)


mafore <- t(ave.foreres$modAveFore)
colnam <- faultname

colnames(mafore) <- colnam

mafore1 <- as.mcmc(mafore[,ind1])
mafore2 <- as.mcmc(mafore[,ind2])
mafore3 <- as.mcmc(mafore[,ind3])


postscript('MAonlyforeDens1Final.eps',paper="special",width=4,height=15,
           onefile = TRUE, horizontal = FALSE,fonts=c("serif", "Palatino"))
par(mar=c(4,5,0.5,0.5))
color_scheme_set("brightblue")
mcmc_areas(mafore1,prob = 0.5, prob_outer = 0.95,point_est = c("median"))
dev.off()


postscript('MAonlyforeDens2Final.eps',paper="special",width=4,height=15,
           onefile = TRUE, horizontal = FALSE,fonts=c("serif", "Palatino"))
par(mar=c(4,5,0.5,0.5))
color_scheme_set("brightblue")
mcmc_areas(mafore2,prob = 0.5, prob_outer = 0.95,point_est = c("median"))
dev.off()

postscript('MAonlyforeDens3Final.eps',paper="special",width=4,height=15,
           onefile = TRUE, horizontal = FALSE,fonts=c("serif", "Palatino"))
par(mar=c(4,5,0.5,0.5))
color_scheme_set("brightblue")
x <- mcmc_areas(log10(mafore3),prob = 0.5, prob_outer = 0.95,point_est = c("median"))
x + scale_x_continuous(breaks = c(3.3,4,5,6,7))
dev.off()
####################################################################







###  Section 3: Clustering or periodic
#########################################################################
###    A regression model for the relationship between the shape 
###    parameter of a Weibull renewal process and the earthquake rate, 
###    faulting type, and tectonic region of each of the 93 faults 
###    suggests that when the earthquake rate increases by 10 fold 
###    (e.g., from 0.0001 to 0.001), the shape parameter of the 
###    Weibull renewal process increases by about 23.6\% (95\% CI: 23.0\%--24.2\%). 
###    The shape parameter of the Weibull renewal process for faults 
###    located at or near plate boundary is about 70% of that for faults 
###    located at an active intraplate, suggesting that the latter appear 
###    to be more periodic than the former. On average, the Weibull 
###    shape parameter for reverse faults is about 31.1\% (95\% CI: 30.2\%-31.9\%)
###    lower than that for normal faults, suggesting that normal faults 
###    appear more periodic than reverse faults.
#########################################################################

##  Earthquake rates

setwd("../ResultsFinal")
## Load the file containing all the parameter estimates
pois.lambda <- weib.alpha <- matrix(0,nfault,15000)
for (faulti in 1:nfault){
  # Poisson process results
  load(paste("B:/Research/PaleoEQ/ResultsFinal/PoisRes/PoisFault-mcmc",faulti,".image",sep=""))

  pois.mat <- as.matrix(mod.mcmc)
  pois.dat <- as.data.frame(pois.mat)
  lambda <- pois.dat$lambda
  pois.lambda[faulti,] <- lambda
  
  # Weibull renewal process results
  load(paste("B:/Research/PaleoEQ/ResultsFinal/WeibRes/WeibFault-mcmc",faulti,".image",sep=""))

  weib.mat <- as.matrix(mod.mcmc)
  weib.dat <- as.data.frame(weib.mat)
  alpha <- weib.dat$alpha
  weib.alpha[faulti,] <- alpha
}  
  
## faulti <- 1
ks.test(log(weib.alpha[faulti,]),"pnorm",mean=mean(log(weib.alpha[faulti,])),sd=sd(log(weib.alpha[faulti,])))
ks.test((weib.alpha[faulti,]),"pnorm",mean=mean((weib.alpha[faulti,])),sd=sd((weib.alpha[faulti,])))

log10rate <- log10(pois.lambda)


meta <- read.csv("../DataFinal/metadata_summary.csv")
attach(meta)
names(meta)
Tec_region <- meta$Tectonic_region
Fault_style <- meta$Faulting_style

Tec_region[Tec_region=="Subduction"] <- "S"
Tec_region[Tec_region=="Plate_boundary_master"] <- "B"
Tec_region[Tec_region=="Plate_boundary_network"] <- "B"
Tec_region[Tec_region=="Intraplate_noncratonic"] <- "I"
Tec_region[Tec_region=="Active_intraplate"] <- "IA"
Tec_region[Tec_region=="Near_plate_boundary"] <- "B"
Tec_region[Tec_region=="Intraplate_cratonic"] <- "I"

Fault_style[Fault_style=="Reverse"] <- "R"
Fault_style[Fault_style=="Strike_slip"] <- "SS"
Fault_style[Fault_style=="Normal"] <- "N"

tecreg12fc <- Tec_region
tecreg12fc[Tec_region=="S"] <- "PB"
tecreg12fc[Tec_region=="B"] <- "PB"
tecreg12fc[Tec_region=="I"] <- "IP"
tecreg12fc[Tec_region=="AI"] <- "IP"
tecreg12fc[Tec_region=="NB"] <- "IP"

Tec_reg <- as.numeric(as.factor(Tec_region))
Fault_st <- as.numeric(as.factor(Fault_style))



## c25 <- glm(log(params$weib.alpha[,2])~Noccur+Tec_reg+Fault_st+log10(rate))

##### Bayesian method lognormal
tmp1 <- "
model{
  for (i in 1:nmcmc){
    for (j in 1:nfault){ # jth fault
     weib.alpha[j,i] ~ dlnorm(mu[j,i],tau)
     mu[j,i] <- a[i] + b[i] * Noccur[j] + c[Tec_reg[j],i] + d[Fault_st[j],i] + g[i] * log10rate[j,i]
   }
   a[i] <- mua + epsa[i]
   epsa[i] ~ dnorm(0,1/asig/asig)
   b[i] <- mub + epsb[i]
   epsb[i] ~ dnorm(0,1/bsig/bsig)
   expb[i] <- exp(b[i])
   for (k in 2:4){
     c[k,i] <- muc[k] + epsc[k,i]
     epsc[k,i] ~ dnorm(0,1/csig[k]/csig[k])
   }
   expc2[i] <- exp(c[2,i])
   expc3[i] <- exp(c[3,i])
   expc4[i] <- exp(c[4,i])
#   expc5[i] <- exp(c[5,i])
   c[1,i] <- 0
   for (kk in 2:3){
     d[kk,i] <- mud[kk] + epsd[kk,i]
     epsd[kk,i] ~ dnorm(0,1/dsig[kk]/dsig[kk])
   }
   expd2[i] <- exp(d[2,i])
   expd3[i] <- exp(d[3,i])
   d[1,i] <- 0
   g[i] <- mug + epsg[i]
   epsg[i] ~ dnorm(0,1/gsig/gsig)
   expg[i] <- exp(g[i])
  }
  mua ~ dnorm(0,1.0E-06)
  asig ~ dt(0,0.04,3)T(0,)
  mub ~ dnorm(0,1.0E-06)
  bsig ~ dt(0,0.04,3)T(0,)
  for (k in 2:4){
    muc[k] ~ dnorm(0,1.0E-06)
    csig[k] ~ dt(0,0.04,3)T(0,)
  }
  for (kk in 2:3){
    mud[kk] ~ dnorm(0,1.0E-06)
    dsig[kk] ~ dt(0,0.04,3)T(0,)
  }
  mug ~ dnorm(0,1.0E-06)
  gsig ~ dt(0,0.04,3)T(0,)
  tau ~ dt(0,0.04,3)T(0,)
}
"

cat(tmp1, file = "weibShapetecreg.jags")

library(R2jags)
nmcmc <- 1000

jagsdata <- list("weib.alpha", "log10rate", "nfault","nmcmc","Noccur","Tec_reg","Fault_st")
params <- c("mua","mub","muc[2]","muc[3]","muc[4]","mud[2]","mud[3]",
            "mug","tau","asig","bsig","csig[2]","csig[3]","csig[4]",
            "dsig[2]","dsig[3]","gsig","expb","expc2","expc3",
            "expc4","expd2","expd3","expg")


system.time(
  weibshape <- jags(data=jagsdata, inits=NULL,
                    parameters.to.save = params,n.chain=3,
                    n.iter=25000, n.thin=20, n.burnin=5000,
                    model.file = 'weibShapetecreg.jags')
)

# Convert to an MCMC object
mod.mcmc <- as.mcmc(weibshape)

save(mod.mcmc,file=paste("weibShapetecreg-mcmc4tecreg.image",sep=""))
#setwd("../ResultsFinal")
#load(paste("weibShapetecreg-mcmc.image",sep=""))

plot(mod.mcmc,trace=TRUE)
gelman.diag(mod.mcmc)
summary(mod.mcmc)

mod.mat <- as.matrix(mod.mcmc)
mod.dat <- as.data.frame(mod.mat)


## Quantile of coefficient exp(mub)
bm <- eval(parse(text=paste('mod.dat$"mub"',sep="")))
expb.quant <- quantile(exp(bm),c(0.025,0.05,0.5,0.95,0.975))
expb.quant
#     2.5%        5%       50%       95%     97.5% 
# 0.9769122 0.9770177 0.9775834 0.9781805 0.9782943


## Quantile of coefficient exp(muc[k])
cm2 <- eval(parse(text=paste('mod.dat$"muc[2]"',sep="")))
expb.quant <- quantile(exp(cm2),c(0.025,0.05,0.5,0.95,0.975))
expb.quant
#    2.5%        5%       50%       95%     97.5% 
# 0.8508384 0.8536225 0.8651944 0.8775757 0.8793929

cm3 <- eval(parse(text=paste('mod.dat$"muc[3]"',sep="")))
expb.quant <- quantile(exp(cm3),c(0.025,0.05,0.5,0.95,0.975))
expb.quant
#     2.5%       5%      50%      95%    97.5% 
# 1.418723 1.420287 1.430026 1.440356 1.441862 

cm4 <- eval(parse(text=paste('mod.dat$"muc[4]"',sep="")))
expb.quant <- quantile(exp(cm4),c(0.025,0.05,0.5,0.95,0.975))
expb.quant
#    2.5%       5%      50%      95%    97.5% 
#  1.837110 1.842611 1.869763 1.897658 1.903555 

## Quantile of coefficient exp(mud[kk])
dm2 <- eval(parse(text=paste('mod.dat$"mud[2]"',sep="")))
expb.quant <- quantile(exp(dm2),c(0.025,0.05,0.5,0.95,0.975))
expb.quant
#    2.5%        5%       50%       95%     97.5% 
# 0.6815740 0.6829019 0.6895461 0.6970615 0.6980696  

dm3 <- eval(parse(text=paste('mod.dat$"mud[3]"',sep="")))
expb.quant <- quantile(exp(dm3),c(0.025,0.05,0.5,0.95,0.975))
expb.quant
#    2.5%        5%       50%       95%     97.5% 
# 0.9046456 0.9059006 0.9117563 0.9180427 0.9190989  

## Quantile of coefficient exp(mug)
gm <- eval(parse(text=paste('mod.dat$"mug"',sep="")))
expb.quant <- quantile(exp(gm),c(0.025,0.05,0.5,0.95,0.975))
expb.quant
#    2.5%       5%      50%      95%    97.5% 
# 1.231516 1.232552 1.237808 1.242634 1.243702







#####################################################################
####        FREQUENTIST METHOD
#####################################################################
load(paste("ParamEstNew.image",sep=""))
rate <- params$pois.lambda[,2]

maxweight <- NULL
for (i in 1:nfault){
  maxweight[i] <- which.max(weights[i,2:6])
}


bmforequant <- Forep50quantile$p50Pois.quant
for (i in 1:nfault){
  if (maxweight[i]==1) bmforequant[i,] <- Forep50quantile$p50Pois.quant[i,]
  if (maxweight[i]==2) bmforequant[i,] <- Forep50quantile$p50gam.quant[i,]
  if (maxweight[i]==3) bmforequant[i,] <- Forep50quantile$p50weib.quant[i,]
  if (maxweight[i]==4) bmforequant[i,] <- Forep50quantile$p50bpt.quant[i,]
  if (maxweight[i]==5) bmforequant[i,] <- Forep50quantile$p50lnorm.quant[i,]
}

bestmod <- NULL
for (i in 1:nfault){
  if (maxweight[i]==1) bestmod[i] <- "P"
  if (maxweight[i]==2) bestmod[i] <- "G"
  if (maxweight[i]==3) bestmod[i] <- "W"
  if (maxweight[i]==4) bestmod[i] <- "B"
  if (maxweight[i]==5) bestmod[i] <- "L"
}


## Log normal regression
Tec_reg <- as.factor(Tec_region)
Fault_st <- as.factor(Fault_style)
mod.bt <- as.factor(bestmod)
c2 <- glm(log(params$weib.alpha[,2])~log10(rate)*mod.bt)
summary(c2)
c21 <- glm(log(params$weib.alpha[,2])~log10(rate))
summary(c22)
c22 <- glm(log(params$weib.alpha[,2])~log10(rate)+mod.bt)
summary(c22)
c23 <- glm(log(params$weib.alpha[,2])~log10(rate)+Noccur)
summary(c23)
c24 <- glm(log(params$weib.alpha[,2])~Noccur)
summary(c24)


##########################################################################
########   SELECTED model

c25 <- glm(log(params$weib.alpha[,2])~Noccur+Tec_reg+Fault_st+log10(rate))
summary(c25)

########   SELECTED model
##########################################################################


AIC(c2);AIC(c21);AIC(c22);AIC(c23);AIC(c24);AIC(c25)


rsd5 <- rstandard(c25)
lp5<-predict(c25)
ks.test(rsd5,"pnorm")


est <- exp(coef(c25))
cf95 <- confint(c25,level=0.95)
cf90 <- confint(c25,level=0.9)
cfex95 <- exp(cf95)
cfex90 <- exp(cf90)
cbind(est,cfex95,cfex90)


                 est     2.5 %     97.5 %       5 %       95 %
(Intercept) 6.2369157 2.9425957 13.2193210 3.3203299 11.7154373
Noccur      0.9775891 0.9538877  1.0018795 0.9576591  0.9979339
Tec_regB    0.7002916 0.5166432  0.9492204 0.5425341  0.9039216
Tec_regI    0.6021636 0.3304496  1.0972961 0.3639187  0.9963794
Tec_regNB   0.6972853 0.5109879  0.9515035 0.5371741  0.9051197
Tec_regS    1.3041242 0.6742584  2.5223862 0.7496987  2.2685646
Fault_stR   0.6899208 0.4575043  1.0404070 0.4887400  0.9739140
Fault_stSS  0.9100786 0.6715208  1.2333840 0.7051552  1.1745542
log10(rate) 1.2365220 1.0398655  1.4703695 1.0692304  1.4299879

#####################################################################
####        FREQUENTIST METHOD
#####################################################################







######################################################################
#######  Figure 4 Weibull shape versus log10rate
######################################################################
Tec_legend <- Tec_region
Tec_legend[Tec_region=="S"] <- "Subduction"
Tec_legend[Tec_region=="B"] <- "(Near) Plate boundary"
Tec_legend[Tec_region=="IA"] <- "Active intraplate"
Tec_legend[Tec_region=="I"] <- "Intraplate"


postscript("WeibShapeParamLogVSlograte.eps",paper="special",width=5*1.6,height=5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,1))
par(mar=c(4.5,4.5,1,1))
plot(log10(rate),log(params$weib.alpha[,2]),xlab="Rate of earthquake occurrence",
     ylab="Weibull shape parameter",cex.axis=1.7,cex.lab=1.7,
       pch=Fault_st,col=Tec_reg,cex=1.5,ylim=c(-1.2,2.1),xlim=c(-6,-1.5),axes=F)
axis(1,(-6):(-2),c(expression(10^-6),expression(10^-5),
                   expression(10^-4),expression(10^-3),
                   expression(10^-2)),cex.axis=1.7)
axis(2,(-1):2,round(exp((-1):2),digits=1),cex.axis=1.7)
box()
legend("bottomright",legend=unique(Faulting_style),pch=unique(Fault_st),cex=1.3,bty="n")
legend("topleft",legend=unique(Tec_legend),col=unique(Tec_reg),pch=16,cex=1.3,bty="n")
dev.off()
######################################################################
#######  Figure 4 Weibull shape versus log10rate
######################################################################








##########################################################################
###  For each fault, we removed the last event in the record to carry 
###  out retrospective forecasts using both the model-averaging and 
###  single-best model approaches. The single-best model remained the 
###  same as that for the full dataset for 54\% of the fault segments 
### (50 out of the 93). It appears that the single-best model is more 
###  likely to change for fault segments with fewer recorded earthquakes.
###  On average, the number of earthquakes along the fault segments 
###  for which the single-best model changed is about 25\% (95\% CI: 6\%-40\%) 
### fewer than that for those fault segments for which the single-best model 
### remained the same. 
###########################################################################
setwd("B:/Research/PaleoEQ/RprogramsNatGeo/SimRetroForecast2")
load(paste("RetroModAveForeRes.image",sep=""))

cut <- 0.95

ind <- 1:nfault
ind.Pois <- ind[ave.foreres$weights[,1]>=cut]  
ind.gam <- ind[ave.foreres$weights[,2]>=cut]  
ind.weib <- ind[ave.foreres$weights[,3]>=cut]  
ind.bpt <- ind[ave.foreres$weights[,4]>=cut]  
ind.lnorm <- ind[ave.foreres$weights[,5]>=cut]  
ind.ma <- ind[-sort(c(ind.Pois,ind.gam,ind.weib,ind.bpt,ind.lnorm))]  

length(ind.Pois)
length(ind.gam)
length(ind.weib)
length(ind.bpt)
length(ind.lnorm)

weights.org <- read.csv("../ResultsFinal/WAICweightsNew.csv")

bestM.retro <- bestM.org <- NULL
for (i in 1:nfault){
  bestM.retro[i] <- which.max(ave.foreres$weights[i,1:5])
  bestM.org[i] <- which.max(weights.org[i,2:6])
}

sum(bestM.retro==bestM.org)
Noccur[bestM.retro==bestM.org]
Noccur[bestM.retro!=bestM.org]

bestM.same <- rep(1,nfault)
bestM.same[bestM.retro!=bestM.org] <- 2


########################################################################
#########  Relation between number of events and whether best model changed
########################################################################

##### Bayesian method
tmp <- "
model{
  for (j in 1:nfault){ # jth fault
    Noccur[j] ~ dnegbin(p[j],size)
    p[j] <- size/(size+mu[j])
    log(mu[j]) <- a + b[bestM.same[j]]
  }
  a ~ dnorm(0,1.0E-06)
  b[1] <- 0
  b[2] ~ dnorm(0,1.0E-06)
  size ~ dunif(0.001,1000)
  theta <- pow(1/mean(p),2)
  scaleparam <- mean((1-p)/p) 
}
"

cat(tmp, file = "NoccurRetroBestMsame.jags")

library(R2jags)
jagsdata <- list("Noccur", "bestM.same", "nfault")
params <- c("a","b[2]","size","theta","scaleparam")

system.time(
  NB.Noccur <- jags(data=jagsdata, inits=NULL,
                    parameters.to.save = params,n.chain=3,
                    n.iter=10000, n.thin=10, n.burnin=5000,
                    model.file = 'NoccurRetroBestMsame.jags')
)

# Convert to an MCMC object
mod.mcmc <- as.mcmc(NB.Noccur)
plot(mod.mcmc,trace=TRUE)
gelman.diag(mod.mcmc)
summary(mod.mcmc)

mod.mat <- as.matrix(mod.mcmc)
mod.dat <- as.data.frame(mod.mat)

## Quantile of coefficient exp(b)
expb.quant <- matrix(NA,nrow=1,ncol=5)
bm <- eval(parse(text=paste('mod.dat$"b[',2,']"',sep="")))
expb.quant <- quantile(exp(bm),c(0.025,0.05,0.5,0.95,0.975))
1-expb.quant







##########################################################################
##  Figure 5 shows the 95\% credible intervals of the forecast of the last 
##  earthquake occurrence time, with 0 representing the mean of the recorded 
##  last earthquake occurrence time. Out of the 93 fault segments
##########################################################################

Fore.Pois <- ave.foreres$Fore.Pois
Fore.gam <- ave.foreres$Fore.gam
Fore.weib <- ave.foreres$Fore.weib
Fore.bpt <- ave.foreres$Fore.bpt
Fore.lnorm <- ave.foreres$Fore.lnorm
modAveFore <- ave.foreres$modAveFore


maxweight <- NULL
for (i in 1:nfault){
  maxweight[i] <- which.max(ave.foreres$weights[i,1:5])
}


bestfore <- t(ave.foreres$modAveFore)
for (i in 1:nfault){
  if (maxweight[i]==1) bestfore[,i] <- (ave.foreres$Fore.Pois)[i,]
  if (maxweight[i]==2) bestfore[,i] <- (ave.foreres$Fore.gam)[i,]
  if (maxweight[i]==3) bestfore[,i] <- (ave.foreres$Fore.weib)[i,]
  if (maxweight[i]==4) bestfore[,i] <- (ave.foreres$Fore.bpt)[i,]
  if (maxweight[i]==5) bestfore[,i] <- (ave.foreres$Fore.lnorm)[i,]
}


colnam <- NULL
for (faulti in 1:nfault)
  colnam <- c(colnam,paste("BM.fore[",faulti,"]",sep=""))

colnames(bestfore) <- colnam


Pois.quant <- gam.quant <- weib.quant <- bpt.quant <- lnorm.quant <- modAv.quant <- bm.quant <- lastocc <- matrix(NA,nrow=nfault,ncol=3)
lastocctimes <- matrix(NA,nrow=nfault,ncol=K)
Pois95width <- MA95width <- NULL
for (faulti in 1:nfault){
  Pois.quant[faulti,] <- quantile(Fore.Pois[faulti,],c(0.025,0.5,0.975))
  gam.quant[faulti,] <- quantile(Fore.gam[faulti,],c(0.025,0.5,0.975))
  weib.quant[faulti,] <- quantile(Fore.weib[faulti,],c(0.025,0.5,0.975))
  bpt.quant[faulti,] <- quantile(Fore.bpt[faulti,],c(0.025,0.5,0.975))
  lnorm.quant[faulti,] <- quantile(Fore.lnorm[faulti,],c(0.025,0.5,0.975))
  modAv.quant[faulti,] <- quantile(modAveFore[faulti,],c(0.025,0.5,0.975))
  bm.quant[faulti,] <- quantile(bestfore[,faulti],c(0.025,0.5,0.975))
  
  Pois95width[faulti] <- Pois.quant[faulti,3]-Pois.quant[faulti,1]
  MA95width[faulti] <- modAv.quant[faulti,3]-modAv.quant[faulti,1]
  
  data = read.csv(paste('../chronologies_all_final/',a[faulti],sep=''), header=FALSE)
  t.occ = data.matrix(data)
  N1 <- ncol(t.occ)
  lastocc[faulti,] <- quantile(t.occ[1:K,N1],c(0.025,0.5,0.975))
  lastocctimes[faulti,] <- t.occ[1:K,N1]
  lastocc[faulti,2] <- mean(lastocctimes[faulti,])
}


x <- modAv.quant - lastocc[,2]
y <- bm.quant - lastocc[,2]

x.p <- Pois.quant - lastocc[,2]
x.g <- gam.quant - lastocc[,2]
x.w <- weib.quant - lastocc[,2]
x.b <- bpt.quant - lastocc[,2]
x.l <- lnorm.quant - lastocc[,2]


sum(x.p[,1]<=0 & x.p[,3]>=0)
sum(x.g[,1]<=0 & x.g[,3]>=0)
sum(x.w[,1]<=0 & x.w[,3]>=0)
sum(x.b[,1]<=0 & x.b[,3]>=0)
sum(x.l[,1]<=0 & x.l[,3]>=0)
sum(x[,1]<=0 & x[,3]>=0)
sum(y[,1]<=0 & y[,3]>=0)

indNotCover <- (1:nfault)[-which(x[,1]<=0 & x[,3]>=0)]
Noccur[indNotCover]


x[indNotCover,]
y[indNotCover,]

(x[ind.ma,3]-x[ind.ma,1])-(y[ind.ma,3]-y[ind.ma,1])


#########################################################################
###  Poisson forecast 95% CI ... times wider than MA forecast 95% CI
#########################################################################
PoisVSma <- Pois95width / MA95width
hist(PoisVSma)
mean(PoisVSma)

## For those fault segments that the MA forecasts covered the true occurrence
mean(PoisVSma[(1:nfault)[-indNotCover]])
hist((PoisVSma[(1:nfault)[-indNotCover]]))

#########################################################################



y1 <- Pois.quant - lastocc[,2]



pow <- 1/4

postscript("RetroMA-PoisCenteredFore1-4th.eps",paper="special",width=8*1.8,height=15*1,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,2))
par(mar=c(4.5,4.1,3.5,0.5))

x.lim <- range(c(range(x),range(y1)))
x.lim <- (abs(x.lim)^(pow))*sign(x.lim)
plot(x.lim,c(0,48),type="n",axes=F,cex.lab=1.8,xlab="Forecast time of last earthquake (year)",
     ylab="",ylim=c(0,48),xlim=c(-20,60))
#axis(1,at=seq(-5,5,1),c(expression(-10^5),expression(-10^4),expression(-10^3),
#                      expression(-10^2),expression(-10),0,10,expression(10^2),
#                      expression(10^3),expression(10^4),expression(10^5)),cex.axis=1.7)
#axis(1,at=c(-5,-3,-1,0,1,3,5),c(expression(-10^5),expression(-10^3),
#                      expression(-10),0,10,expression(10^3),expression(10^5)),cex.axis=1.7)
axis(1,at=c(-(10^5)^(pow),-(10^4)^(pow),-100^(pow),0,100^(pow),(10^4)^(pow),
            (10^6)^(pow),(10^7)^(pow)),
     c(expression(-10^5),expression(-10^4),
       "",0,"",expression(10^4),expression(10^6),
       expression(10^7)),cex.axis=1.7)
axis(3,at=c(-(10^5)^(pow),-(10^4)^(pow),-100^(pow),0,100^(pow),(10^4)^(pow),
            (10^6)^(pow),(10^7)^(pow)),
     c("","",
       expression(-10^2),"",expression(10^2),"","",
       ""),cex.axis=1.7)
axis(2,at=seq(1,47,1),seq(1,47,1),cex.axis=1.7)

for (faulti in 1:47){
  
  temma <- (modAv.quant[faulti,]- lastocc[faulti,2])
  tempois <- (Pois.quant[faulti,] - lastocc[faulti,2])
  
  points(((abs(temma)^(pow))*sign(temma))[2],faulti,cex=1.8,pch=16,col=1)
  arrows(x0=((abs(temma)^(pow))*sign(temma))[1], y0=faulti, 
         x1=((abs(temma)^(pow))*sign(temma))[3], y1=faulti, code=3, angle=90, 
         length=0.07, col="black", lwd=2)
  
  points(((abs(tempois)^(pow))*sign(tempois))[2],faulti-0.3,cex=1.8,pch=16,col="grey")
  arrows(x0=((abs(tempois)^(pow))*sign(tempois))[1], y0=faulti-0.3, 
         x1=((abs(tempois)^(pow))*sign(tempois))[3], y1=faulti-0.3, code=3, 
         angle=90, length=0.07, col="grey", lwd=2)
  
  #  lines((abs(temma)^(pow))*sign(temma),rep(faulti,3),lwd=1.8,col=1)
  #  lines((abs(tempois)^(pow))*sign(tempois),rep(faulti-0.2,3),lwd=1.8,col="grey")
} 
abline(v=0,col="red",lwd=1.5)


par(mar=c(4.5,1,3.5,0.5))

plot(x.lim,c(0,48),type="n",axes=F,cex.lab=1.8,xlab="Forecast time of last earthquake",
     ylab="",ylim=c(0,48),xlim=c(-20,60))
axis(1,at=c(-(10^5)^(pow),-(10^4)^(pow),-100^(pow),0,100^(pow),(10^4)^(pow),
            (10^6)^(pow),(10^7)^(pow)),
     c(expression(-10^5),expression(-10^4),
       "",0,"",expression(10^4),expression(10^6),
       expression(10^7)),cex.axis=1.7)
axis(3,at=c(-(10^5)^(pow),-(10^4)^(pow),-100^(pow),0,100^(pow),(10^4)^(pow),
            (10^6)^(pow),(10^7)^(pow)),
     c("","",
       expression(-10^2),"",expression(10^2),"","",
       ""),cex.axis=1.7)
axis(2,at=seq(1,46,1),seq(48,93,1),cex.axis=1.7)

for (faulti in 48:nfault){
  
  temma <- (modAv.quant[faulti,]- lastocc[faulti,2])
  tempois <- (Pois.quant[faulti,] - lastocc[faulti,2])
  
  points(((abs(temma)^(pow))*sign(temma))[2],faulti-47,cex=1.8,pch=16,col=1)
  arrows(x0=((abs(temma)^(pow))*sign(temma))[1], y0=faulti-47, 
         x1=((abs(temma)^(pow))*sign(temma))[3], y1=faulti-47, code=3, angle=90, 
         length=0.07, col="black", lwd=2)
  
  points(((abs(tempois)^(pow))*sign(tempois))[2],faulti-47-0.3,cex=1.8,pch=16,col="grey")
  arrows(x0=((abs(tempois)^(pow))*sign(tempois))[1], y0=faulti-47-0.3, 
         x1=((abs(tempois)^(pow))*sign(tempois))[3], y1=faulti-47-0.3, code=3, 
         angle=90, length=0.07, col="grey", lwd=2)
  
  #  lines((abs(temma)^(pow))*sign(temma),rep(faulti-47,3),lwd=1.8,col=1)
  #  lines((abs(tempois)^(pow))*sign(tempois),rep(faulti-47-0.2,3),lwd=1.8,col="grey")
} 
abline(v=0,col="red",lwd=1.5)
dev.off()























#####################################################################
###   Supplement Figures 1-3 forecasts from all models (since 2022)
#####################################################################
load(paste("B:/Research/PaleoEQ/RprogramsNatGeo/ModAveForeResNew.image",sep=""))

Fore.Pois <- ave.foreres$Fore.Pois
Fore.gam <- ave.foreres$Fore.gam
Fore.weib <- ave.foreres$Fore.weib
Fore.bpt <- ave.foreres$Fore.bpt
Fore.lnorm <- ave.foreres$Fore.lnorm
modAveFore <- ave.foreres$modAveFore


postscript("waicAvFore1.eps",paper="special",width=8,height=12*1.05,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(11,3))
par(mar=c(2,1.5,2.5,1.5))

#for (faulti in ((1:93)[-indnores])[1:47]){
for (faulti in 1:31){
  Pois.quant <- quantile(Fore.Pois[faulti,],c(0.025,0.5,0.975))
  gam.quant <- quantile(Fore.gam[faulti,],c(0.025,0.5,0.975))
  weib.quant <- quantile(Fore.weib[faulti,],c(0.025,0.5,0.975))
  bpt.quant <- quantile(Fore.bpt[faulti,],c(0.025,0.5,0.975))
  lnorm.quant <- quantile(Fore.lnorm[faulti,],c(0.025,0.5,0.975))
  modAv.quant <- quantile(modAveFore[faulti,],c(0.025,0.5,0.975))
  
  x.lim <- range(c(Pois.quant,gam.quant,weib.quant,bpt.quant,lnorm.quant,modAv.quant))
  
  plot(x.lim,c(0,1),type="n",axes=F,cex.lab=1.5,xlab="Forecast time of next earthquake (year)",
       ylab="",main=faultname[faulti],ylim=c(-0.1,1.1))
  axis(1,cex.axis=1.5)
  box()
  points(Pois.quant[2],1,cex=1.8,pch=16,col=6)
  arrows(x0=Pois.quant[1], y0=1, 
         x1=Pois.quant[3], y1=1, code=3, angle=90, 
         length=0.07, col=6, lwd=2)
  
  points(gam.quant[2],0.8,cex=1.8,pch=16,col=2)
  arrows(x0=gam.quant[1], y0=0.8, 
         x1=gam.quant[3], y1=0.8, code=3, angle=90, 
         length=0.07, col=2, lwd=2)

  points(weib.quant[2],0.6,cex=1.8,pch=16,col=3)
  arrows(x0=weib.quant[1], y0=0.6, 
         x1=weib.quant[3], y1=0.6, code=3, angle=90, 
         length=0.07, col=3, lwd=2)

  points(bpt.quant[2],0.4,cex=1.8,pch=16,col=4)
  arrows(x0=bpt.quant[1], y0=0.4, 
         x1=bpt.quant[3], y1=0.4, code=3, angle=90, 
         length=0.07, col=4, lwd=2)
  
  points(lnorm.quant[2],0.2,cex=1.8,pch=16,col=5)
  arrows(x0=lnorm.quant[1], y0=0.2, 
         x1=lnorm.quant[3], y1=0.2, code=3, angle=90, 
         length=0.07, col=5, lwd=2)

  points(modAv.quant[2],0,cex=1.8,pch=16,col=1)
  arrows(x0=modAv.quant[1], y0=0, 
         x1=modAv.quant[3], y1=0, code=3, angle=90, 
         length=0.07, col=1, lwd=2)
} 

plot(c(2028.902,3967.561),c(0,1),type="n",axes=F,xlab="",
     ylab="",ylim=c(-0.1,1.1))

xleg <- c(rep(2000,3),rep(2800,3))
yleg <- c(1.2,0.8,0.4,1.2,0.8,0.4)+0.1
leg <- c("Poisson","Gamma","Weibull","BPT","lognormal","MA")
pchleg = c(16,16,16,16,16,16)
colleg = c(6,2,3,4,5,1)
for (i in 1:6)
  legend(xleg[i],yleg[i],leg[i],pch=pchleg[i],col=colleg[i],bty="n",cex=1.5)

dev.off()



postscript("waicAvFore2.eps",paper="special",width=8,height=12*1.05,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(11,3))
par(mar=c(2,1.5,2.5,1.5))

#for (faulti in ((1:93)[-indnores])[1:47]){
for (faulti in 32:62){
  Pois.quant <- quantile(Fore.Pois[faulti,],c(0.025,0.5,0.975))
  gam.quant <- quantile(Fore.gam[faulti,],c(0.025,0.5,0.975))
  weib.quant <- quantile(Fore.weib[faulti,],c(0.025,0.5,0.975))
  bpt.quant <- quantile(Fore.bpt[faulti,],c(0.025,0.5,0.975))
  lnorm.quant <- quantile(Fore.lnorm[faulti,],c(0.025,0.5,0.975))
  modAv.quant <- quantile(modAveFore[faulti,],c(0.025,0.5,0.975))
  
  x.lim <- range(c(Pois.quant,gam.quant,weib.quant,bpt.quant,lnorm.quant,modAv.quant))
  
  plot(x.lim,c(0,1),type="n",axes=F,cex.lab=1.5,xlab="Forecast time of next earthquake (year)",
       ylab="",main=faultname[faulti],ylim=c(-0.1,1.1))
  axis(1,cex.axis=1.5)
  box()
  points(Pois.quant[2],1,cex=1.8,pch=16,col=6)
  arrows(x0=Pois.quant[1], y0=1, 
         x1=Pois.quant[3], y1=1, code=3, angle=90, 
         length=0.07, col=6, lwd=2)
  
  points(gam.quant[2],0.8,cex=1.8,pch=16,col=2)
  arrows(x0=gam.quant[1], y0=0.8, 
         x1=gam.quant[3], y1=0.8, code=3, angle=90, 
         length=0.07, col=2, lwd=2)
  
  points(weib.quant[2],0.6,cex=1.8,pch=16,col=3)
  arrows(x0=weib.quant[1], y0=0.6, 
         x1=weib.quant[3], y1=0.6, code=3, angle=90, 
         length=0.07, col=3, lwd=2)
  
  points(bpt.quant[2],0.4,cex=1.8,pch=16,col=4)
  arrows(x0=bpt.quant[1], y0=0.4, 
         x1=bpt.quant[3], y1=0.4, code=3, angle=90, 
         length=0.07, col=4, lwd=2)
  
  points(lnorm.quant[2],0.2,cex=1.8,pch=16,col=5)
  arrows(x0=lnorm.quant[1], y0=0.2, 
         x1=lnorm.quant[3], y1=0.2, code=3, angle=90, 
         length=0.07, col=5, lwd=2)
  
  points(modAv.quant[2],0,cex=1.8,pch=16,col=1)
  arrows(x0=modAv.quant[1], y0=0, 
         x1=modAv.quant[3], y1=0, code=3, angle=90, 
         length=0.07, col=1, lwd=2)
} 

plot(c(2028.902,3967.561),c(0,1),type="n",axes=F,xlab="",
     ylab="",ylim=c(-0.1,1.1))

xleg <- c(rep(2000,3),rep(2800,3))
yleg <- c(1.2,0.8,0.4,1.2,0.8,0.4)+0.1
leg <- c("Poisson","Gamma","Weibull","BPT","lognormal","MA")
pchleg = c(16,16,16,16,16,16)
colleg = c(6,2,3,4,5,1)
for (i in 1:6)
  legend(xleg[i],yleg[i],leg[i],pch=pchleg[i],col=colleg[i],bty="n",cex=1.5)

dev.off()


postscript("waicAvFore3.eps",paper="special",width=8,height=12*1.05,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(11,3))
par(mar=c(2,1.5,2.5,1.5))

#for (faulti in ((1:93)[-indnores])[1:47]){
for (faulti in 63:93){
  Pois.quant <- quantile(Fore.Pois[faulti,],c(0.025,0.5,0.975))
  gam.quant <- quantile(Fore.gam[faulti,],c(0.025,0.5,0.975))
  weib.quant <- quantile(Fore.weib[faulti,],c(0.025,0.5,0.975))
  bpt.quant <- quantile(Fore.bpt[faulti,],c(0.025,0.5,0.975))
  lnorm.quant <- quantile(Fore.lnorm[faulti,],c(0.025,0.5,0.975))
  modAv.quant <- quantile(modAveFore[faulti,],c(0.025,0.5,0.975))
  
  x.lim <- range(c(Pois.quant,gam.quant,weib.quant,bpt.quant,lnorm.quant,modAv.quant))
  
  plot(x.lim,c(0,1),type="n",axes=F,cex.lab=1.5,xlab="Forecast time of next earthquake (year)",
       ylab="",main=faultname[faulti],ylim=c(-0.1,1.1))
  axis(1,cex.axis=1.5)
  box()
  points(Pois.quant[2],1,cex=1.8,pch=16,col=6)
  arrows(x0=Pois.quant[1], y0=1, 
         x1=Pois.quant[3], y1=1, code=3, angle=90, 
         length=0.07, col=6, lwd=2)
  
  points(gam.quant[2],0.8,cex=1.8,pch=16,col=2)
  arrows(x0=gam.quant[1], y0=0.8, 
         x1=gam.quant[3], y1=0.8, code=3, angle=90, 
         length=0.07, col=2, lwd=2)
  
  points(weib.quant[2],0.6,cex=1.8,pch=16,col=3)
  arrows(x0=weib.quant[1], y0=0.6, 
         x1=weib.quant[3], y1=0.6, code=3, angle=90, 
         length=0.07, col=3, lwd=2)
  
  points(bpt.quant[2],0.4,cex=1.8,pch=16,col=4)
  arrows(x0=bpt.quant[1], y0=0.4, 
         x1=bpt.quant[3], y1=0.4, code=3, angle=90, 
         length=0.07, col=4, lwd=2)
  
  points(lnorm.quant[2],0.2,cex=1.8,pch=16,col=5)
  arrows(x0=lnorm.quant[1], y0=0.2, 
         x1=lnorm.quant[3], y1=0.2, code=3, angle=90, 
         length=0.07, col=5, lwd=2)
  
  points(modAv.quant[2],0,cex=1.8,pch=16,col=1)
  arrows(x0=modAv.quant[1], y0=0, 
         x1=modAv.quant[3], y1=0, code=3, angle=90, 
         length=0.07, col=1, lwd=2)
} 

plot(c(2028.902,3967.561),c(0,1),type="n",axes=F,xlab="",
     ylab="",ylim=c(-0.1,1.1))

xleg <- c(rep(2000,3),rep(2800,3))
yleg <- c(1.2,0.8,0.4,1.2,0.8,0.4)+0.1
leg <- c("Poisson","Gamma","Weibull","BPT","lognormal","MA")
pchleg = c(16,16,16,16,16,16)
colleg = c(6,2,3,4,5,1)
for (i in 1:6)
  legend(xleg[i],yleg[i],leg[i],pch=pchleg[i],col=colleg[i],bty="n",cex=1.5)

dev.off()












#####################################################################
###   Supplement Figure MSE Retro forecast
#####################################################################

mseres <- read.csv("B:/Research/PaleoEQ/RprogramsNatGeo/SimRetroForecast2/MAforeBiasVarMSE.csv")

MAbias <- mseres$MAbias
poisbias <- mseres$poisbias
gambias <- mseres$gambias
weibbias <- mseres$weibbias
bptbias <- mseres$bptbias
lnormbias <- mseres$lnormbias

MAvar <- mseres$MAvar
poisvar <- mseres$poisvar
gamvar <- mseres$gamvar
weibvar <- mseres$weibvar
bptvar <- mseres$bptvar
lnormvar <- mseres$lnormvar

MAmse <- mseres$MAmse
poismse <- mseres$poismse
gammse <- mseres$gammse
weibmse <- mseres$weibmse
bptmse <- mseres$bptmse
lnormmse <- mseres$lnormmse

rangebias <- range(c(MAbias,poisbias,gambias,weibbias,bptbias,lnormbias))
rangevar <- range(c(MAvar,poisvar,gamvar,weibvar,bptvar,lnormvar))
rangemse <- range(c(MAmse,poismse,gammse,weibmse,bptmse,lnormmse))


#####################################################################
#############################  NOT USED   #################################
postscript("RetroMA-BMforeMSE1-10th.eps",paper="special",width=5*cos (35.4/180*pi)/0.612,height=5*3,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(3,1))
par(mar=c(4.5,4.5,4.5,1))
plot(1:93,(abs(MAbias)^(1/10))*sign(MAbias),ylim=c((abs(rangebias[1])^(1/10))*sign(rangebias[1]),
                                                   (rangebias[2])^(1/10)),cex=1.5,cex.axis=1.5,cex.lab=1.5,xlab="Fault number (MA)",
     ylab="Bias",pch=16)
points(1:93,(abs(poisbias)^(1/10))*sign(poisbias),cex=1.5,pch=15,col=2)
#axis(3,1:93,Noccur,cex.axis=1.5)

par(mar=c(4.5,4.5,4.5,1))
plot(1:93,(MAvar)^(1/10),ylim=c((rangevar[1])^(1/10),20),
     cex=1.5,cex.axis=1.5,cex.lab=1.5,xlab="Fault number (MA)",
     ylab="Variance",pch=16)
points(1:93,(poisvar)^(1/10),cex=1.5,pch=15,col=2)
#axis(3,1:93,Noccur,cex.axis=1.5)

par(mar=c(4.5,4.5,4.5,1))
plot(1:93,(MAmse)^(1/10),ylim=c(0,20),cex=1.5,cex.axis=1.5,cex.lab=1.5,xlab="Fault number (MA)",
     ylab="Mean squared error",pch=16)
points(1:93,(poismse)^(1/10),cex=1.5,pch=15,col=2)
#axis(3,1:93,Noccur,cex.axis=1.5)
xleg <- c(40,80)
yleg <- rep(18,2)
leg <- c("MA","Poisson")
pchleg = c(16,15)
colleg = 1:2
for (i in 1:2)
  legend(xleg[i],yleg[i],leg[i],pch=pchleg[i],col=colleg[i],bty="n",cex=1.8)
dev.off()
#############################  NOT USED   #################################
#####################################################################


MSEdiff <- log10(abs(poismse-MAmse))*sign(poismse-MAmse)
MSEdiff[is.na(MSEdiff)] <- 0
postscript("RetroMA-BMforeMSEVSnoccur.eps",paper="special",
           width=5*cos (35.4/180*pi)/0.612,height=5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,1))
par(mar=c(4.5,4.5,1,1))
#hist(log10(abs(poismse-MAmse))*sign(poismse-MAmse),main="",
#     xlab="Difference in mean squared error",ylab="Frequency",cex.axis=1.7,
#     cex.lab=1.7,axes=F)
#axis(1,at=c(-10,-5,0,5,10),c(expression(-10^10),expression(-10^5),
#                             0,expression(10^5),expression(10^10)),cex.axis=1.7)

plot(Noccur-1,MSEdiff,xlab="Number of earthquakes",
     ylab="Difference in MSE (Poisson-MA)",cex.axis=1.7,cex.lab=1.7,
     axes=F)
axis(1,cex.axis=1.7)
axis(2,at=c(-10,-5,0,5,10),c(expression(-10^10),expression(-10^5),
                             0,expression(10^5),expression(10^10)),cex.axis=1.7)
box()
dev.off()






 
  
  