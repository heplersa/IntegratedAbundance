
library(coda) # For manipulation of MCMC results.
library(scoringRules)
############################################################
###########             Load data               ############
############################################################


### time intervals
TT=10
n=100
H=100

### load in data
load("SimulatedData.Rda")
Ntrue = N[1:(n*TT),]
ftrue = f[1:(n*TT),]
gtrue = g[1:(n*TT),]
utrue = u[1:(n*TT),]
vtrue = v[1:(n*TT),]
lambdatrue = exp(utrue+vtrue)

filename1 = "JointYearly."
filename2 = "JointSparse."
filename3 = "SingleYearly."
filename4 = "SingleSparse."
filenamelist = c(filename1,filename2,filename3,filename4)

CPall = matrix(0,H,4)
MSEall = matrix(0,H,4)
MSE0all = matrix(0,H,4)
rate.MSEall = matrix(0,H,4)
rate.MSE0all = matrix(0,H,4)
rel.MAEall = matrix(0,H,4)
rel.MAE0all = matrix(0,H,4)
CRPSall = matrix(0,H,4)
CRPS0all = matrix(0,H,4)
beta0all = matrix(0,H,4)
beta1all = matrix(0,H,4)

#output all intercept estimates for just the first case
betaThat0 = matrix(0,TT,H)
betaDhat0 = matrix(0,TT,H)


for(i in 1:4){
  
  Nhat=matrix(0,n*TT,H)
  lambdahat = matrix(0,n*TT,H)
  Nsd=matrix(0,n*TT,H)
  uhat = matrix(0,n*TT,H)
  fhat=matrix(0,n*TT,H)
  ghat = matrix(0,n*TT,H)
  Nlb = matrix(0,n*TT,H)
  Nub = matrix(0,n*TT,H)
  betahat = matrix(0,2,H)
  gammahat = matrix(0,2,H)
  betaThat = matrix(0,TT,H)
  betaDhat = matrix(0,TT,H)
  MSE = rep(0,H)  ### overall MSE
  MSE.cont = matrix(0,n*TT,H)
  rate.MSE = rep(0,H)
  rate.MSE0 = rep(0,H)
  rel.MAE = rep(0,H)
  MSE0 = rep(0,H)
  rel.MAE0 = rep(0,H)
  CP = rep(0,H)
  S.hat = matrix(0,TT,H)
  CRPS = rep(0,H)
  CRPS0 = rep(0,H)

  
  if(i %in% c(2,4)){
    for(h in 1:H){
      for(t in 1:TT){
        temp = lm(S[c(2,5,8),h]~c(2,5,8))
        S.hat[t,h] = temp$coefficients[1]+t*temp$coefficients[2]
      }
    }
  }else{
    for(h in 1:H){
      for(t in 1:TT){
        temp = lm(S[,h]~c(1:TT))
        S.hat[t,h] = temp$coefficients[1]+t*temp$coefficients[2]
      }
    }
  }
  
  for(m in c(1:4)){
    filesname = paste(filenamelist[i], toString(m), ".Rda",sep="")
    load(filesname)
    results=SimMCMC
    kk=0
    for(k in ((m-1)*25+1):(25*m)){
      kk=kk+1
      ll = min(which(names(results[[kk]][[1]])=="N[1]"))
      uu = min(which(names(results[[kk]][[1]])=="N[1000]"))
      Nhat[,k] = results[[kk]][[1]][ll:uu]
      Nlb[,k] = results[[kk]][[2]][1,ll:uu]
      Nub[,k] = results[[kk]][[2]][2,ll:uu]
      ll = min(which(names(results[[kk]][[1]])=="lambda[1]"))
      uu = min(which(names(results[[kk]][[1]])=="lambda[1000]"))
      lambdahat[,k] = results[[kk]][[1]][ll:uu]
      ll = min(which(names(results[[kk]][[1]])=="beta[1]"))
      uu = min(which(names(results[[kk]][[1]])=="beta[2]"))
      betahat[,k] = results[[kk]][[1]][ll:uu]
      ll = min(which(names(results[[kk]][[1]])=="beta.D[1]"))
      uu = min(which(names(results[[kk]][[1]])=="beta.D[10]"))
      betaDhat[,k] = results[[kk]][[1]][ll:uu]
      if(i<3){
        ll = min(which(names(results[[kk]][[1]])=="beta.T[1]"))
        uu = min(which(names(results[[kk]][[1]])=="beta.T[10]"))
        betaThat[,k] = results[[kk]][[1]][ll:uu]
      }
      MSE[k] = sqrt(mean((Nhat[,k]-Ntrue[,k])^2))
      MSE.cont[,k] = (Nhat[,k]-Ntrue[,k])^2
      MSE0[k] = sqrt(mean((kronecker(S.hat[,k],rep(1,n))*P-Ntrue[,k])^2))
      rate.MSE[k] = sqrt(mean((lambdahat[,k]-lambdatrue[,k])^2))
      rate.MSE0[k] = sqrt(mean((rep(1,1000)-lambdatrue[,k])^2))
      rel.MAE[k] = median(abs(Nhat[,k]-Ntrue[,k])/Ntrue[,k])
      rel.MAE0[k] = median(abs(kronecker(S.hat[,k],rep(1,n))*P-Ntrue[,k])/Ntrue[,k])
      CP[k] = mean(Ntrue[,k]>Nlb[,k] & Ntrue[,k]<Nub[,k])
      beta0all[k,i]=betahat[1,k]
      beta1all[k,i]=betahat[2,k]
      CRPS[k] = mean(crps_binom(Ntrue[,k],P, kronecker((betahat[1,k]+betahat[2,k]*(1:TT)),rep(1,n))*lambdahat[,k]))
      CRPS0[k] = mean(crps_binom(Ntrue[,k],P,kronecker(S.hat[,k],rep(1,n))))
    }
  }
  
  if(i==1){
    betaThat0=betaThat
    betaDhat0=betaDhat
  }
  CPall[,i] = CP
  MSEall[,i] = MSE
  MSE0all[,i] = MSE0
  rate.MSEall[,i] = rate.MSE
  rate.MSE0all[,i] = rate.MSE0
  rel.MAEall[,i] = rel.MAE
  rel.MAE0all[,i] = rel.MAE0
  CRPSall[,i]=CRPS
  CRPS0all[,i] = CRPS0
}


###### boxplot with proposed all, proposed sparse, baseline all
png('ProposedModelBoxplots.png',width=16,height=16,res=300,units='in')
par(mar=c(6,4.1,4.1,2.1))
par(mgp=c(3,1.5,0))
par(mfrow=c(2,2))
boxplot(cbind(MSEall[,1],MSEall[,2],MSE0all[,1]),main = "Count RMSE",names=c("Yearly", "Sparse","Baseline"),cex.axis=3,pch=16,cex=1.2,cex.main=3,cex.lab=3)
boxplot(cbind(rate.MSEall[,1],rate.MSEall[,2],rate.MSE0all[,1]),main="Relative Rate RMSE", names=c("Yearly", "Sparse","Baseline"),cex.axis=3,pch=16,cex=1.2,cex.main=3,cex.lab=3)
boxplot(cbind(rel.MAEall[,1],rel.MAEall[,2],rel.MAE0all[,1]),main = "Count Relative MAE",names=c("Yearly", "Sparse","Baseline"),cex.axis=3,pch=16,cex=1.2,cex.main=3,cex.lab=3)
boxplot(cbind(CRPSall[,1],CRPSall[,2],CRPS0all[,1]),main = "Count RPS",names=c("Yearly", "Sparse","Baseline"),cex.axis=3,pch=16,cex=1.2,cex.main=3,cex.lab=3)
dev.off()


png('ProposedandDeath.png',width=16,height=16,res=300,units='in')
par(mar=c(6,4.1,6.1,2.1))
par(mgp=c(3,1.5,0))
par(cex.axis=3)
par(cex.main=3)
par(cex.lab=3)
par(mfrow=c(2,2))
boxplot(cbind(MSEall[,1],MSEall[,2],MSEall[,3],MSEall[,4]),main= "Count RMSE", names=c("Yearly","Sparse","Yearly","Sparse"),pch=16,cex=1.2)
xtick = c(1.5,3.5)
xticklab=c("Joint Model","Single Outcome")
axis(side=1, at=xtick, tick=FALSE, labels = xticklab,padj=1.75)
boxplot(cbind(rate.MSEall[,1],rate.MSEall[,2],rate.MSEall[,3],rate.MSEall[,4]),main= "Relative Rate RMSE", names=c("Yearly","Sparse","Yearly","Sparse"),pch=16,cex=1.2)
axis(side=1, at=xtick,  tick=FALSE,labels = xticklab,padj=1.75)
boxplot(cbind(rel.MAEall[,1],rel.MAEall[,2],rel.MAEall[,3],rel.MAEall[,4]),main= "Count Relative MAE", names=c("Yearly","Sparse","Yearly","Sparse"),pch=16,cex=1.2)
axis(side=1, at=xtick,  tick=FALSE,labels = xticklab,padj=1.75)
boxplot(cbind(CRPSall[,1],CRPSall[,2],CRPSall[,3],CRPSall[,4]),main= "Count RPS", names=c("Yearly","Sparse","Yearly","Sparse"),pch=16,cex=1.2)
axis(side=1, at=xtick,  tick=FALSE,labels = xticklab,padj=1.75)
dev.off()


####### plots of intercept parameters

hist(beta0all[,1],main="",xlab=expression(beta[0]^{mu}))
abline(v=0.05,col="red",lwd=2)

hist(beta1all[,1],main="",xlab=expression(beta[1]^{mu}))
abline(v=-0.0001,col="red",lwd=2)


par(mar=c(5.1, 4.8, 4.1, 2.1))
par(mgp=c(3,1,0))
boxplot(t(betaThat0),xlab="Year",ylab=expression(mu[t]^{(1)}))
lines(seq(0.5,4.5,.5),rep(-3,length(seq(0.5,4.5,.5))),col="red",lwd=2)
lines(seq(4.5,7.5,.5),rep(-2.5,length(seq(4.5,7.5,.5))),col="red",lwd=2)
lines(seq(7.5,10.5,.5),rep(-2,length(seq(7.5,10.5,.5))),col="red",lwd=2)

boxplot(t(betaDhat0),xlab="Year",ylab=expression(mu[t]^{(2)}))
lines(seq(0.5,5.5,.5),rep(-5,length(seq(0.5,5.5,.5))),col="red",lwd=2)
lines(seq(5.5,7.5,.5),rep(-4.5,length(seq(5.5,7.5,.5))),col="red",lwd=2)
lines(seq(7.5,10.5,.5),rep(-4,length(seq(7.5,10.5,.5))),col="red",lwd=2)


################################################################################
###### Stacked boxplots for all simulated data sets 

filename = "JointSparse."

Nhat=matrix(0,n*TT,H)
lambdahat = matrix(0,n*TT,H)

for(h in 1:H){
  for(t in 1:TT){
    temp = lm(S[c(2,5,8),h]~c(2,5,8))
    S.hat[t,h] = temp$coefficients[1]+t*temp$coefficients[2]
  }
}


for(m in c(1:4)){
  filesname = paste(filename, toString(m), ".Rda",sep="")
  load(filesname)
  results=SimMCMC
  kk=0
  for(k in ((m-1)*25+1):(25*m)){
    kk=kk+1
    ll = min(which(names(results[[kk]][[1]])=="N[1]"))
    uu = min(which(names(results[[kk]][[1]])=="N[1000]"))
    Nhat[,k] = results[[kk]][[1]][ll:uu]
    ll = min(which(names(results[[kk]][[1]])=="lambda[1]"))
    uu = min(which(names(results[[kk]][[1]])=="lambda[1000]"))
    lambdahat[,k] = results[[kk]][[1]][ll:uu]
  }
}

### random choose 3 data sets to make scatterplots
set.seed(1132021)
Hind = sort(sample(1:H,3,replace=FALSE))

ctplot = list()
rtplot = list()
htplot = list()

for(i in 1:3){
  datanew = data.frame("Nhat"=Nhat[,Hind[i]],"Pop"=P,"Ntrue"=Ntrue[,Hind[i]],"lambdahat"=lambdahat[,Hind[i]],"lambdatrue"=lambdatrue[,Hind[i]])
  ctplot[[i]] = ggplot(datanew,aes(x=Nhat,y=Ntrue)) +
    geom_point() +
    geom_abline(intercept=0,slope=1,color="red",size=1.1)+
    labs(x="Estimated Count", y="True Count", title=paste("Simulation ", Hind[i]))+
    theme(text=element_text(size=16))
  
  rtplot[[i]] = ggplot(datanew,aes(x=lambdahat,y=lambdatrue)) +
    geom_point() +
    geom_abline(intercept=0,slope=1,color="red",size=1.1)+
    labs(x="Estimated Rate", y="True Rate", title=paste("Simulation ", Hind[i]))+
    theme(text=element_text(size=16))
  
  htplot[[i]] = ggplot(datanew,aes(lambdatrue-lambdahat))+
    geom_histogram() + theme(text=element_text(size=16))+
    labs(x="Estimated Difference Between\nTrue Rate and Estimated Rate", title=paste("Simulation ", Hind[i]))
  
}

#library(gridExtra)

g=grid.arrange(arrangeGrob(ctplot[[1]],ctplot[[2]],ctplot[[3]],nrow=1),widths=unit(11,"in"),heights=unit(4,"in"))
ggsave("CountScatterplot.png",plot=g,units="in") 

g2=grid.arrange(arrangeGrob(rtplot[[1]],rtplot[[2]],rtplot[[3]],nrow=1),widths=unit(11,"in"),heights=unit(4,"in"))
ggsave("RateScatterplot.png",plot=g2,units="in") 

gh=grid.arrange(arrangeGrob(htplot[[1]],htplot[[2]],htplot[[3]],nrow=1),widths=unit(11,"in"),heights=unit(4,"in"))
ggsave("RateDiffHistogram.png",plot=gh,units="in") 

#### now make side by side boxplots of the difference in rates for all simulations with vertical lines for quartiles

Nhatdata = data.frame("Simulation" = as.factor(kronecker(1:H,rep(1,1000))),"Estimates"=as.vector(Nhat),"Population"=as.vector(P),"Actual"=as.vector(Ntrue),"ActualRate"=as.vector(lambdatrue),"EstimatedRate"=as.vector(lambdahat))
Nhatdata$diff = Nhatdata$ActualRate-Nhatdata$EstimatedRate
Nhatdata$diff2 = Nhatdata$Actual/Nhatdata$Population-Nhatdata$Estimates/Nhatdata$Population
Nhatdata$relbias = (Nhatdata$ActualRate-Nhatdata$EstimatedRate)/Nhatdata$ActualRate

gp=ggplot(data=Nhatdata,aes(x=Simulation,y=(diff)))+
  geom_boxplot(outlier.size=.2,outlier.alpha=.5,fill="gray")+coord_flip()+
  scale_y_continuous(name="Difference Between True and Estimated Relative Risk",breaks=seq(-2,4.25,by=0.5))+
  geom_hline(aes(yintercept=quantile(diff,c(0.025))),color="red",linetype="solid")+
  geom_hline(aes(yintercept=quantile(diff,c(0.5))),color="red",linetype="solid")+
  geom_hline(aes(yintercept=quantile(diff,c(0.975))),color="red",linetype="solid")+
  scale_x_discrete(labels=NULL,breaks=NULL)+ theme(text=element_text(size=16))

ggsave("RateDiffBoxplots.png",plot=gp,units="in") 

gp2=ggplot(data=Nhatdata,aes(x=Simulation,y=relbias))+
  geom_boxplot(outlier.size=.5,outlier.alpha=.5,fill="gray")+coord_flip()+
  scale_y_continuous(name="Relative Difference Betwen True and Estimated Relative Risk",breaks=seq(-3,1,by=0.5))+
  geom_hline(aes(yintercept=quantile(relbias,c(0.025))),color="red",linetype="solid")+
  geom_hline(aes(yintercept=quantile(relbias,c(0.5))),color="red",linetype="solid")+
  geom_hline(aes(yintercept=quantile(relbias,c(0.975))),color="red",linetype="solid")+
  scale_x_discrete(labels=NULL,breaks=NULL)+ theme(text=element_text(size=16))

ggsave("RelBiasBoxplots.png",plot=gp2,units="in") 


