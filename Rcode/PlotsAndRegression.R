
###################################################################
###########   Code for all plots and regression results  ##########
###################################################################

library(rjags)
a <- dir("../DataFinal/chronologies_all_final")
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



Noccur <- NULL  ## a vector of number of earthquakes along each fault segment
mean.inter <- sd.inter <- NULL
X.div.Xave <- list()
sd.X.div.Xave <- matrix(NA,nrow=nfault,ncol=K)
sd.X.div.Xave.med <- matrix(NA,nrow=nfault,ncol=3)

#######################################################################
#####  Scaled inter-event times: inter-event times X/X_ave        #####
#####     X_ave is the mean inter-event times) for each fault     #####
#######################################################################
for (faulti in (1:nfault)){
  data = read.csv(paste('../DataFinal/chronologies_all_final/',a[faulti],sep=''), header=FALSE)

# Chronologies have occurrence times in calendar years. 
# m: number of MC samples of the chronologies
# N: Number of earthquakes at the fault
  t.occ = data.matrix(data)
  m = nrow(t.occ)
  N <- ncol(t.occ)
  
## Interarrival times mean and variance
  diffoccur <- apply(t.occ[1:K,],1,diff)
  mean.inter <- append(mean.inter,mean(diffoccur))
  sd.inter <- append(sd.inter,sd(diffoccur))

  diff.mean <- apply(diffoccur,2,mean)
  diff.sd <- apply(diffoccur,2,sd)
  
  temp <- diffoccur
  for (kk in 1:K) temp[,kk] <- diffoccur[,kk]/diff.mean[kk]
  X.div.Xave[[faulti]] <- temp
  
  sd.X.div.Xave[faulti,] <- apply(X.div.Xave[[faulti]],2,sd)

  sd.X.div.Xave.med[faulti,] <- quantile(sd.X.div.Xave[faulti,],c(0.025,0.5,0.975))

  Noccur[faulti] <- N  
}



###############################################################
#####################   Read Results   ########################
###############################################################
load(paste("../Results/ModAveProbRes.image",sep=""))
load(paste("../Results/ModAveLogProbRes.image",sep=""))
load(paste("../Results/ModAveForeRes.image",sep=""))

weights <- read.csv("../Results/WAICweights.csv")

## Load the file containing the quantiles of the parameter estimates
load(paste("../Results/ParamEst.image",sep=""))


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


## Save the quantiles of each forecast of probability of an earthquake in the next 50 years on each fault.
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





#####################################################################
###  Supplementary Fig. 4 SD of X/Xave for faults with the same bestM
#####################################################################

postscript("XdivXaveSDbestM-rate.eps",paper="special",width=5*1.6*1.2,height=5*1.5*1.2,
           onefile = TRUE, horizontal = FALSE)
par(mfrow=c(3,2))
par(mar=c(4.5,4.5,1,1))

sd.X.div.Xave.weib <- sd.X.div.Xave.med[ind.weib,]
sd.X.div.Xave.weib.sort <- sd.X.div.Xave.weib[order(sd.X.div.Xave.weib[,2]),]
plot(1:length(ind.weib),sd.X.div.Xave.weib.sort[,1],ylim=c(0,1.3),type="n",
     xlab="Index",
     ylab="SD of scaled inter-event time",cex.axis=1.7,cex.lab=1.7)
for (i in 1:length(ind.weib)){
  points(rep(i,3),sd.X.div.Xave.weib.sort[i,],col=c(2,1,2))
}
matlines(1:length(ind.weib),sd.X.div.Xave.weib.sort,col=c(2,1,2),lty=c(2,1,2))
legend("bottomright",legend=c("Median","95% CI"),col=c(1,2),lty=c(1,2),cex=1.5,bty="n",pch=1)
legend("topleft",legend="Weibull",cex=1.5,bty="n")
text(20,1.2,"(a)",cex=1.5)

sd.X.div.Xave.gam <- sd.X.div.Xave.med[ind.gam,]
sd.X.div.Xave.gam.sort <- sd.X.div.Xave.gam[order(sd.X.div.Xave.gam[,2]),]
plot(1:length(ind.gam),sd.X.div.Xave.gam.sort[,1],ylim=c(0,1.3),type="n",
     xlab="Index",
     ylab="SD of scaled inter-event time",cex.axis=1.7,cex.lab=1.7)
for (i in 1:length(ind.gam)){
  points(rep(i,3),sd.X.div.Xave.gam.sort[i,],col=c(2,1,2))
}
matlines(1:length(ind.gam),sd.X.div.Xave.gam.sort,col=c(2,1,2),lty=c(2,1,2))
legend("bottomright",legend=c("Median","95% CI"),col=c(1,2),lty=c(1,2),cex=1.5,bty="n",pch=1)
legend("topleft",legend="Gamma",cex=1.5,bty="n")
text(3,1.2,"(b)",cex=1.5)


sd.X.div.Xave.bpt <- sd.X.div.Xave.med[ind.bpt,]
sd.X.div.Xave.bpt.sort <- sd.X.div.Xave.bpt[order(sd.X.div.Xave.bpt[,2]),]
plot(1:length(ind.bpt),sd.X.div.Xave.bpt.sort[,1],ylim=c(0,2.5),type="n",
     xlab="Index",
     ylab="SD of scaled inter-event time",cex.axis=1.7,cex.lab=1.7)
for (i in 1:length(ind.bpt)){
  points(rep(i,3),sd.X.div.Xave.bpt.sort[i,],col=c(2,1,2))
}
matlines(1:length(ind.bpt),sd.X.div.Xave.bpt.sort,col=c(2,1,2),lty=c(2,1,2))
legend("bottomright",legend=c("Median","95% CI"),col=c(1,2),lty=c(1,2),cex=1.5,bty="n",pch=1)
legend("topleft",legend="BPT",cex=1.5,bty="n")
text(2.5,2.4,"(c)",cex=1.5)


sd.X.div.Xave.lnorm <- sd.X.div.Xave.med[ind.lnorm,]
sd.X.div.Xave.lnorm.sort <- sd.X.div.Xave.lnorm[order(sd.X.div.Xave.lnorm[,2]),]

plot(1:length(ind.lnorm),sd.X.div.Xave.lnorm.sort[,1],ylim=c(0,2.1),type="n",
     xlab="Index",
     ylab="SD of scaled inter-event time",cex.axis=1.7,cex.lab=1.7)
for (i in 1:length(ind.lnorm)){
  points(rep(i,3),sd.X.div.Xave.lnorm.sort[i,],col=c(2,1,2))
}
matlines(1:length(ind.lnorm),sd.X.div.Xave.lnorm.sort,col=c(2,1,2),lty=c(2,1,2))
legend("bottomright",legend=c("Median","95% CI"),col=c(1,2),lty=c(1,2),cex=1.5,bty="n",pch=1)
legend("topleft",legend="lognormal",cex=1.5,bty="n")
text(8,2,"(d)",cex=1.5)



pois.lambda <- params$pois.lambda[,2]
log10rate <- log10(pois.lambda)

sd.X.div.median <- sd.X.div.Xave.med[,2]

plot(log10rate,(sd.X.div.median),xlab="Rate of earthquake occurrence",
     ylab="SD of scaled inter-event time",cex.axis=1.7,cex.lab=1.7,
     cex=1.5,xlim=c(-6,-1.5),axes=F)
axis(1,(-6):(-2),c(expression(10^-6),expression(10^-5),
                   expression(10^-4),expression(10^-3),
                   expression(10^-2)),cex.axis=1.7)
axis(2,(-1):2,round(exp((-1):2),digits=1),cex.axis=1.7)
text(-3.5,log(9),"(e)",cex=1.5)
box()
dev.off()




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
###   best model have on average about 1.59 (90% CI (1.10,2.47)) times 
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

########################################################################
#########  Relation between number of events and best model
########################################################################

##### Bayesian method
NoccurBestMod <- "
model{
  for (j in 1:nfault){ # jth fault
    Noccur[j] ~ dpois(mu[j])
    log(mu[j]) <- a + b[X[j]]
  }
  a ~ dnorm(0,1.0E-06)
  b[1] <- 0
  for (i in 2:5){
    b[i] ~ dnorm(0,1.0E-06)
  }
}
"


library(R2jags)
X <- mod.best.no
jagsdata <- list("Noccur", "X", "nfault")
params <- c("a","b[2]","b[3]","b[4]","b[5]")

system.time(
  BM.Noccur <- jags(data=jagsdata, inits=NULL,
                    parameters.to.save = params,n.chain=3,
                    n.iter=30000, n.thin=20, n.burnin=10000,
                    model.file = textConnection(NoccurBestMod))
)

# Convert to an MCMC object
mod.mcmc <- as.mcmc(BM.Noccur)
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

#b2 Gamm  [1,] 1.0858542 1.1633290 1.701171 2.597468 2.837502
#b3 Lnorm [2,] 0.7668302 0.8019029 1.140877 1.688493 1.835715
#b4 MA    [3,] 0.7718876 0.8156574 1.125864 1.651279 1.776576
#b5 Weib  [4,] 1.1032292 1.1646982 1.590568 2.293427 2.471833
#########################################################################





# Section 2

########################################################################
###############    Table 1
####  Print fault name and MA forecast and Best model forecast
########################################################################

table1 <- cbind(faultname,rep("&",93),Noccur,rep("&",93),
     formatC(100*Forep50quantile$p50modAve.quant[,2], format = "f", digits = 3),
     rep("(",93),
     formatC(100*Forep50quantile$p50modAve.quant[,1], format = "f", digits = 3),
     rep(",",93),
     formatC(100*Forep50quantile$p50modAve.quant[,3], format = "f", digits = 3),
     rep(")",93),rep("&",93),
     bestmod,rep("&",93),
     formatC(as.matrix(100*bmforequant[,2]), format = "f", digits = 3),
     rep("(",93),
     formatC(as.matrix(100*bmforequant[,1]), format = "f", digits = 3),
     rep(",",93),
     formatC(as.matrix(100*bmforequant[,3]), format = "f", digits = 3),rep(")",93))

write.table(table1,"Table1.txt")
write.csv(table1,"Table1.csv")




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
######   On average, faults along a plate boundary are about 41 (95% CI 
######   (10,162)) times more likely to have a large earthquake in 
######   the next 50 years than intraplate faults.
#########################################################################
###  Load meta data on tectonic settings and faulting style.
meta <- read.csv("../DataFinal/metadata_summary.csv")
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
tecreg12 <- as.numeric(as.factor(tecreg12fc))
unique(tecreg12fc)
unique(tecreg12)


### Load probability of an event in the next 50 years
# load("../Results/ModAveProbRes.image")
p50modAve.mean <- p50modAve.sd <- NULL
lp50modAve.mean <- lp50modAve.sd <- NULL
for (faulti in 1:nfault){
  p50modAve.vec <- as.vector(ave.probres$p50modAve[faulti,,])+1e-8
  p50modAve.mean[faulti] <- mean(p50modAve.vec)
  p50modAve.sd[faulti] <- sd(p50modAve.vec)
  lp50modAve.sd[faulti] <- sd(log(p50modAve.vec))
}

##### Bayesian method truncated lognormal (0,1)
p50lnorm <- "
model{
    for (j in 1:nfault){ # jth fault
      # Model as follows
      # 1. Estimated value ~ logNormal( True Value (log), Estimation error (log scale))
      p50modAve.mean[j] ~ dlnorm(lp50[j], 1/lp50modAve.sd[j]/lp50modAve.sd[j])T(0,1)
      
      # 2. (log) True value ~ Normal(Process average, Process variation) 
      lp50[j] ~ dnorm(mu[j], tau.process)
      
      # 3. Deterministic Process average 
      mu[j] <- a + b[tecreg12[j]]
   }
   a ~ dnorm(0,1.0E-06)
   b[2] ~ dnorm(0,1.0E-06)
   expb <- exp(b[2])
   b[1] <- 0

   tau.process ~ dt(0,0.04,3)T(0,)
   sd.process <- 1/sqrt(tau.process)
}
"

library(R2jags)
jagsdata <- list("p50modAve.mean", "nfault","tecreg12","lp50modAve.sd")
params <- c("a","b[2]","expb", "sd.process")

system.time(
  p50.tec <- jags(data=jagsdata, inits=NULL,
                   parameters.to.save = params,n.chain=3,
                   n.iter=55000, n.thin=20, n.burnin=5000,
                   model.file = textConnection(p50lnorm))
)

# Convert to an MCMC object
mod.mcmc <- as.mcmc(lp50.tec)


plot(mod.mcmc,trace=TRUE)
gelman.diag(mod.mcmc)
summary(mod.mcmc)

mod.mat <- as.matrix(mod.mcmc)
mod.dat <- as.data.frame(mod.mat)


## Quantile of coefficient exp(mub)

## Quantile of coefficient exp(muc[k])
bm2 <- eval(parse(text=paste('mod.dat$"b[2]"',sep="")))
expb.quant <- quantile(exp(bm2),c(0.025,0.05,0.5,0.7,0.95,0.975))
expb.quant 
## 2.5%        5%       50%       70%       95%     97.5% 
##10.38133  12.91447  41.06738  58.87504 129.97939 162.17024 

expb <- eval(parse(text=paste('mod.dat$expb',sep="")))
expb.quant <- quantile(expb,c(0.025,0.05,0.5,0.7,0.95,0.975))
expb.quant 

#####################################################################








########################################################################
## Difference between averaged and best model 50 year forecast probability
########################################################################
which(abs(Forep50quantile$p50modAve.quant[,3]-bmforequant[,3])>=0.001)
#[1] 31 35 42 61 67 76 81 91
which(abs(Forep50quantile$p50modAve.quant[,1]-bmforequant[,1])>=0.001)
#[1] 11 12 14 30 44 61 64 67 76 81

sort(abs(Forep50quantile$p50modAve.quant[,3]-bmforequant[,3]))
sort(abs(Forep50quantile$p50modAve.quant[,2]-bmforequant[,2]))
sort(abs(Forep50quantile$p50modAve.quant[,1]-bmforequant[,1]))

sum(abs(Forep50quantile$p50modAve.quant[,3]-bmforequant[,3])<=0.001)/93
sum(abs(Forep50quantile$p50modAve.quant[,2]-bmforequant[,2])<=0.001)/93
sum(abs(Forep50quantile$p50modAve.quant[,1]-bmforequant[,1])<=0.001)/93
########################################################################

########################################################################
## Difference between averaged and best model forecast next event times
## See after Figure 3
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





########################################################################
## Difference between averaged and best model forecast next-event times
########################################################################
ma.fore.quant <- t(apply(t(ave.foreres$modAveFore),2,quantile,probs =c(0.025,0.5,0.975)))

maxweight <- NULL
for (i in 1:nfault){
  maxweight[i] <- which.max(weights[i,2:6])
}

bestfore <- t(ave.foreres$modAveFore)
for (i in 1:nfault){
  if (maxweight[i]==1) bestfore[,i] <- t(ave.foreres$Fore.Pois)[,i]
  if (maxweight[i]==2) bestfore[,i] <- t(ave.foreres$Fore.gam)[,i]
  if (maxweight[i]==3) bestfore[,i] <- t(ave.foreres$Fore.weib)[,i]
  if (maxweight[i]==4) bestfore[,i] <- t(ave.foreres$Fore.bpt)[,i]
  if (maxweight[i]==5) bestfore[,i] <- t(ave.foreres$Fore.lnorm)[,i]
}

bm.fore.quant <- t(apply(bestfore,2,quantile,probs =c(0.025,0.5,0.975)))


sort(abs(ma.fore.quant[,3]-bm.fore.quant[,3]))
sort(abs(ma.fore.quant[,2]-bm.fore.quant[,2]))
sort(abs(ma.fore.quant[,1]-bm.fore.quant[,1]))

max(abs(ma.fore.quant[,3]-bm.fore.quant[,3]))
max(abs(ma.fore.quant[,2]-bm.fore.quant[,2]))
max(abs(ma.fore.quant[,1]-bm.fore.quant[,1]))

sum(abs(ma.fore.quant[,3]-bm.fore.quant[,3])[ind.ma]>=50)
sum(abs(ma.fore.quant[,2]-bm.fore.quant[,2])[ind.ma]>=50)
sum(abs(ma.fore.quant[,1]-bm.fore.quant[,1])[ind.ma]>=50)
########################################################################









###  Section 3: Clustering or periodic


######################################################################
###  Based on the parameter estimates from the Gamma and Weibull 
###  renewal processes, five fault segments appear to show clustering 
###  behaviour, with both of the upper 95% credible limit of the shape 
###  parameter from each of these models being less than 1. These are 
###  Cadell, Dunstan, Lake Edgar, Solitario Canyon, and Waitangi, all 
###  of which have low earthquake occurrence rates. For Dunstan and 
###  Waitangi, the BPT renewal process has a WAIC weight over 0.95, 
###  and the estimates of $\beta$ are over 2.5, confirming a clustering 
###  behaviour. Six fault segments appear to demonstrate near Poisson 
###  behaviour, Dead Sea Beteiha, Dead Sea Qatar, Dead Sea Taybeh, 
###  Langshan Piedmont Xibulong East, Reelfoot, and Wharekuri, 
###  with the 95\% credible interval of the shape parameter for 
###  both the Gamma and Weibull model containing 1.
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





#########################################################################
###    A regression model for the relationship between the shape parameter 
###    of a Weibull renewal process and the earthquake rate, tectonic 
###    setting, faulting type, and the number of earthquakes of each of 
###.   the 93 fault segments suggests that when the earthquake rate increases 
###.   10 fold (e.g., from 0.0001 to 0.001), the shape parameter of the 
###.   Weibull renewal process increases by about 23.8% (95% CI: 5.0%--45.4%).
###.   ...On average, the Weibull shape parameter for reverse faults is 
###    about 69.0% (95% CI: 45.9%-104.3%) of that for normal faults, 
###.   suggesting that large earthquakes recur more periodically on normal 
###.   faults than on reverse faults. 
#########################################################################

##  Earthquake rates
## Load the file containing all the parameter estimates
pois.lambda <- weib.alpha <- matrix(0,nfault,15000)
pois.lambda.mean <- weib.alpha.mean <- weib.alpha.sd <- weib.lalpha.mean <- NULL
weib.lalpha.sd <- weib.lalpha.smeansd <- NULL
for (faulti in 1:nfault){
  # Poisson process results
  load(paste("../Results/PoisRes/PoisFault-mcmc",faulti,".image",sep=""))
  
  pois.mat <- as.matrix(mod.mcmc)
  pois.dat <- as.data.frame(pois.mat)
  lambda <- pois.dat$lambda
  pois.lambda[faulti,] <- lambda
  pois.lambda.mean[faulti] <- mean(lambda)
  
  # Weibull renewal process results
  load(paste("../Results/WeibRes/WeibFault-mcmc",faulti,".image",sep=""))
  
  weib.mat <- as.matrix(mod.mcmc)
  weib.dat <- as.data.frame(weib.mat)
  alpha <- weib.dat$alpha
  weib.alpha[faulti,] <- alpha
  weib.alpha.mean[faulti] <- mean(alpha)
  weib.alpha.sd[faulti] <- sd(alpha)
  weib.lalpha.sd[faulti] <- sd(log(alpha))
  weib.lalpha.mean[faulti] <- mean(log(alpha))
  weib.lalpha.smeansd[faulti] <- sd(log(alpha))/sqrt(length(alpha))
}  


log10rate <- log10(pois.lambda.mean)


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

Tec_reg <- as.numeric(as.factor(Tec_region))
Fault_st <- as.numeric(as.factor(Fault_style))



##### Bayesian method log alpha ~ normal
weibShapetecreg <- "
model{
   for (j in 1:nfault){ # jth fault
      # Model as follows
      # 1. (log) Estimated value ~ Normal((log) True Value, Estimation error (log scale))
      weib.lalpha.mean[j] ~ dnorm(weib.lalpha[j], 1/weib.lalpha.sd[j]/weib.lalpha.sd[j])
      
      # 2. (log) True value ~ Normal(Process average, Process variation) 
      weib.lalpha[j] ~ dnorm(mu[j], tau.process)
      
      # 3. Deterministic Process average 
      mu[j] <- a + b * Noccur[j] + c[Tec_reg[j]] + d[Fault_st[j]] + g * log10rate[j]
   }
   a ~ dnorm(0,1.0E-06)
   b ~ dnorm(0,1.0E-06)
   expb <- exp(b)
   for (k in 2:4){
     c[k] ~ dnorm(0,1.0E-06)
     expc[k] <- exp(c[k])
   }
   c[1] <- 0
   expc[1] <- 1
 
   for (kk in 2:3){
     d[kk] ~ dnorm(0,1.0E-06)
     expd[kk] <- exp(d[kk])
   }
   d[1] <- 0
   expd[1] <- 1
   g ~ dnorm(0,1.0E-06)
   expg <- exp(g)
 
   tau.process ~ dt(0,0.04,3)T(0,)
   sd.process <- 1/sqrt(tau.process)
}
"


library(R2jags)

jagsdata <- list("weib.lalpha.mean", "log10rate", "nfault","Noccur","Tec_reg",
                 "Fault_st","weib.lalpha.sd")
params <- c("a","b","c[2]","c[3]","c[4]","d[2]","d[3]",
            "g","expb","expc[2]","expc[3]","expc[4]",#"tau",
            "expd[2]","expd[3]","expg", "sd.process")


system.time(
  weibshape <- jags(data=jagsdata, inits=NULL,
                    parameters.to.save = params,n.chain=3,
                    n.iter=55000, n.thin=20, n.burnin=5000,
                    model.file = textConnection(weibShapetecreg))
)

# Convert to an MCMC object
mod.mcmc <- as.mcmc(weibshape)


plot(mod.mcmc,trace=TRUE)
gelman.diag(mod.mcmc)
summary(mod.mcmc)

mod.mat <- as.matrix(mod.mcmc)
mod.dat <- as.data.frame(mod.mat)


## Quantile of coefficient exp(mub)
bm <- eval(parse(text=paste('mod.dat$"b"',sep="")))
expb.quant <- quantile(exp(bm),c(0.025,0.05,0.5,0.95,0.975))
expb.quant
#2.5%        5%       50%       95%     97.5% 
#0.9540273 0.9578144 0.9776678 0.9980505 1.0021928 

## Quantile of coefficient exp(muc[k])
cm2 <- eval(parse(text=paste('mod.dat$"c[2]"',sep="")))
expb.quant <- quantile(exp(cm2),c(0.025,0.05,0.5,0.7,0.95,0.975))
expb.quant 
#2.5%        5%       50%       95%     97.5% 
#0.5061693 0.5509567 0.8631802 1.3511186 1.4755947 

cm3 <- eval(parse(text=paste('mod.dat$"c[3]"',sep="")))
expb.quant <- quantile(exp(cm3),c(0.025,0.05,0.5,0.95,0.975))
expb.quant
#2.5%       5%      50%      95%    97.5% 
#1.098021 1.148646 1.434963 1.789391 1.867041 

cm4 <- eval(parse(text=paste('mod.dat$"c[4]"',sep="")))
expb.quant <- quantile(exp(cm4),c(0.025,0.05,0.5,0.95,0.975))
expb.quant
#2.5%       5%      50%      95%    97.5% 
#1.043492 1.158048 1.857226 3.023668 3.318071 

## Quantile of coefficient exp(mud[kk])
dm2 <- eval(parse(text=paste('mod.dat$"d[2]"',sep="")))
expb.quant <- quantile(exp(dm2),c(0.025,0.05,0.5,0.95,0.975))
expb.quant
#2.5%        5%       50%       95%     97.5% 
#0.4593140 0.4914297 0.6902697 0.9783300 1.0435317 

dm3 <- eval(parse(text=paste('mod.dat$"d[3]"',sep="")))
expb.quant <- quantile(exp(dm3),c(0.025,0.05,0.5,0.95,0.975))
expb.quant
#2.5%        5%       50%       95%     97.5% 
#0.6922591 0.7242296 0.9120356 1.1450206 1.1984509 

## Quantile of coefficient exp(mug)
gm <- eval(parse(text=paste('mod.dat$"g"',sep="")))
expb.quant <- quantile(exp(gm),c(0.025,0.05,0.5,0.95,0.975))
expb.quant
#2.5%       5%      50%      95%    97.5% 
#1.050303 1.077166 1.238058 1.414628 1.454355 





#####################################################################
###  Supplementary Fig. 5 Posterior distributions of the parameters of the regression mode
#####################################################################

postscript("PosteriorDensWeibShapeReg.eps",paper="special",width=5*1.6*1.5,height=5*1.6,
           onefile = TRUE, horizontal = FALSE)
par(mfrow=c(4,2))
par(mar=c(4.5,4.5,1.5,1))
densplot(mod.mcmc[,15],cex.axis=1.7,cex.lab=1.7,
         main="Parameter for log10(rate)",xlab=expression(e^g),ylab="Density",lwd=2)
densplot(mod.mcmc[,10],cex.axis=1.7,cex.lab=1.7,
         main="Parameter for Intraplate setting",xlab=expression(e^c[2]),ylab="Density",lwd=2)
densplot(mod.mcmc[,11],cex.axis=1.7,cex.lab=1.7,
         main="Parameter for Active Intraplate setting",xlab=expression(e^c[3]),ylab="Density",lwd=2)
densplot(mod.mcmc[,12],cex.axis=1.7,cex.lab=1.7,
         main="Parameter for Subduction setting",xlab=expression(e^c[4]),ylab="Density",lwd=2)
densplot(mod.mcmc[,13],cex.axis=1.7,cex.lab=1.7,
         main="Parameter for Reverse fault",xlab=expression(e^d[2]),ylab="Density",lwd=2)
densplot(mod.mcmc[,14],cex.axis=1.7,cex.lab=1.7,
         main="Parameter for Strike slip fault",xlab=expression(e^d[3]),ylab="Density",lwd=2)
densplot(mod.mcmc[,9],cex.axis=1.7,cex.lab=1.7,
         main="Parameter for No of events",xlab=expression(e^b),ylab="Density",lwd=2)
dev.off()






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
plot(log10rate,weib.lalpha.mean,xlab="Rate of earthquake occurrence",
     ylab="Weibull shape parameter",cex.axis=1.7,cex.lab=1.7,
     pch=as.numeric(Fault_st)-1,col=as.numeric(Tec_reg),
     ylim=c(-1.2,2.1),xlim=c(-6,-1.5),axes=F)
axis(1,(-6):(-2),c(expression(10^-6),expression(10^-5),
                   expression(10^-4),expression(10^-3),
                   expression(10^-2)),cex.axis=1.7)
axis(2,(-1):2,round(exp((-1):2),digits=1),cex.axis=1.7)
box()
legend("bottomright",legend=unique(Faulting_style),pch=unique(as.numeric(Fault_st)-1),cex=1.3,bty="n")
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
weights.retro <- read.csv("../Results/RetroWAICweights.csv")

cut <- 0.95

ind <- 1:nfault
ind.Pois <- ind[weights.retro[,2]>=cut]  
ind.gam <- ind[weights.retro[,3]>=cut]  
ind.weib <- ind[weights.retro[,4]>=cut]  
ind.bpt <- ind[weights.retro[,5]>=cut]  
ind.lnorm <- ind[weights.retro[,6]>=cut]  
ind.ma <- ind[-sort(c(ind.Pois,ind.gam,ind.weib,ind.bpt,ind.lnorm))]  

length(ind.Pois)
length(ind.gam)
length(ind.weib)
length(ind.bpt)
length(ind.lnorm)

weights.org <- read.csv("../Results/WAICweights.csv")

bestM.retro <- bestM.org <- NULL
for (i in 1:nfault){
  bestM.retro[i] <- which.max(weights.retro[i,2:6])
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
NoccurRetroBestMsame <- "
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


library(R2jags)
jagsdata <- list("Noccur", "bestM.same", "nfault")
params <- c("a","b[2]","size","theta","scaleparam")

system.time(
  NB.Noccur <- jags(data=jagsdata, inits=NULL,
                    parameters.to.save = params,n.chain=3,
                    n.iter=10000, n.thin=10, n.burnin=5000,
                    model.file = textConnection(NoccurRetroBestMsame))
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
#2.5%        5%       50%       95%     97.5% 
#0.3939125 0.3770268 0.2659697 0.1370057 0.1050767 






##########################################################################
##  Figure 5 shows the 95\% credible intervals of the forecast of the last 
##  earthquake occurrence time, with 0 representing the mean of the recorded 
##  last earthquake occurrence time. Out of the 93 fault segments
##########################################################################
load(paste("../Results/RetroModAveForeRes.image",sep=""))

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
Pois95width <- MA95width <- BM95width <- NULL
for (faulti in 1:nfault){
  Pois.quant[faulti,] <- quantile(Fore.Pois[faulti,],c(0.025,0.5,0.975))
  gam.quant[faulti,] <- quantile(Fore.gam[faulti,],c(0.025,0.5,0.975))
  weib.quant[faulti,] <- quantile(Fore.weib[faulti,],c(0.025,0.5,0.975))
  bpt.quant[faulti,] <- quantile(Fore.bpt[faulti,],c(0.025,0.5,0.975))
  lnorm.quant[faulti,] <- quantile(Fore.lnorm[faulti,],c(0.025,0.5,0.975))
  modAv.quant[faulti,] <- quantile(modAveFore[faulti,],c(0.025,0.5,0.975))
  bm.quant[faulti,] <- quantile(bestfore[,faulti],c(0.025,0.5,0.975))
  
  Pois95width[faulti] <- Pois.quant[faulti,3]-Pois.quant[faulti,1]
  BM95width[faulti] <- bm.quant[faulti,3]-bm.quant[faulti,1]
  MA95width[faulti] <- modAv.quant[faulti,3]-modAv.quant[faulti,1]
  
  data = read.csv(paste('../DataFinal/chronologies_all_final/',a[faulti],sep=''), header=FALSE)
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
Noccur[indNotCover]-1


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
###   Supplementary Fig. 1-3 forecasts from all models (since 2022)
#####################################################################
load(paste("../Results/ModAveForeRes.image",sep=""))

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
###   Supplementary Fig MSE Retro forecast
#####################################################################

mseres <- read.csv("../Results/MAforeBiasVarMSE.csv")

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




MSEdiff <- log10(abs(poismse-MAmse))*sign(poismse-MAmse)
MSEdiff[is.na(MSEdiff)] <- 0

sum(MSEdiff>0)

sum((MSEdiff[Noccur==5]>0))

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





## Relative difference
MSEreldiff <- poismse/MAmse
sum(MSEreldiff>=2)

postscript("RetroMA-PoisforeMSErelativeVSnoccur.eps",paper="special",
           width=5*cos (35.4/180*pi)/0.612,height=5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,1))
par(mar=c(4.5,4.5,1,1))
#hist(log10(abs(poismse-MAmse))*sign(poismse-MAmse),main="",
#     xlab="Difference in mean squared error",ylab="Frequency",cex.axis=1.7,
#     cex.lab=1.7,axes=F)
#axis(1,at=c(-10,-5,0,5,10),c(expression(-10^10),expression(-10^5),
#                             0,expression(10^5),expression(10^10)),cex.axis=1.7)

plot(Noccur-1,log(MSEreldiff),xlab="Number of earthquakes",
     ylab="MSE Poisson/MSE MA",cex.axis=1.7,cex.lab=1.7,ylim=c(log(0.1),log(41)),
     axes=F)
axis(1,cex.axis=1.7)
axis(2,at=c(log(0.1),0,log(10),log(2),log(40)),c(0.1,1,10,2,40),cex.axis=1.7)
box()
dev.off()




maxweight <- NULL
for (i in 1:nfault){
  maxweight[i] <- which.max(ave.foreres$weights[i,1:5])
}

bestbias <- NULL 
bestvar <- NULL 
bestmse <- NULL 
for (i in 1:nfault){
  if (maxweight[i]==1) bestbias[i] <- poisbias[i]
  if (maxweight[i]==2) bestbias[i] <- gambias[i]
  if (maxweight[i]==3) bestbias[i] <- weibbias[i]
  if (maxweight[i]==4) bestbias[i] <- bptbias[i]
  if (maxweight[i]==5) bestbias[i] <- lnormbias[i]
}

for (i in 1:nfault){
  if (maxweight[i]==1) bestvar[i] <- poisvar[i]
  if (maxweight[i]==2) bestvar[i] <- gamvar[i]
  if (maxweight[i]==3) bestvar[i] <- weibvar[i]
  if (maxweight[i]==4) bestvar[i] <- bptvar[i]
  if (maxweight[i]==5) bestvar[i] <- lnormvar[i]
}

for (i in 1:nfault){
  if (maxweight[i]==1) bestmse[i] <- poismse[i]
  if (maxweight[i]==2) bestmse[i] <- gammse[i]
  if (maxweight[i]==3) bestmse[i] <- weibmse[i]
  if (maxweight[i]==4) bestmse[i] <- bptmse[i]
  if (maxweight[i]==5) bestmse[i] <- lnormmse[i]
}



## Relative difference
MSEreldiff <- bestmse/MAmse
postscript("RetroMA-BMforeMSErelativeVSnoccur.eps",paper="special",
           width=5*cos (35.4/180*pi)/0.612,height=5,onefile = TRUE, horizontal = FALSE)
par(mfrow=c(1,1))
par(mar=c(4.5,4.5,1,1))
#hist(log10(abs(poismse-MAmse))*sign(poismse-MAmse),main="",
#     xlab="Difference in mean squared error",ylab="Frequency",cex.axis=1.7,
#     cex.lab=1.7,axes=F)
#axis(1,at=c(-10,-5,0,5,10),c(expression(-10^10),expression(-10^5),
#                             0,expression(10^5),expression(10^10)),cex.axis=1.7)

plot(Noccur-1,log(MSEreldiff),xlab="Number of earthquakes",
     ylab="MSE Best Model/MSE MA",cex.axis=1.7,cex.lab=1.7,ylim=c(log(0.5),log(2)),
     axes=F)
axis(1,cex.axis=1.7)
axis(2,at=c(log(0.5),0,log(10),log(2),log(0.1)),c(0.5,1,10,2,0.1),cex.axis=1.7)
box()
dev.off()

sort((MSEreldiff[Noccur==5]))











  