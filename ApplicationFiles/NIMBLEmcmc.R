
########################################################################################################################
########################## Implement MCMC algorithm for the Ohio data application ######################################
########################################################################################################################

library(nimble) # For MCMC computation using NIMBLE.
library(coda) # For manipulation of MCMC results.

############################################################
###########             Load data               ############
############################################################

data = read.csv(file="ApplicationData.csv",header=T)
T = length(unique(data$Year))
n = length(unique(data$County))

### load the 5 and 1 year ACS data
dataest5=read.csv(file="OH_ACS_5yrEst.csv",header=T)
datamoe5=read.csv(file="OH_ACS_5yrMOE.csv",header=T)
dataest1=read.csv(file="OH_ACS_1yrEst.csv",header=T)
datamoe1=read.csv(file="OH_ACS_1yrMOE.csv",header=T)
### note there is one more county in the 1 year data for 2015-2018 than earlier years. Remove.
#table(dataest1$GEOID) shows FIPS 39009 starts in 2015. Remove this
rmInd = which(dataest1$GEOID=="39009")
dataest1 = dataest1[-rmInd,]
datamoe1 = datamoe1[-rmInd,]

### create an indicator of what counties we have one year data for
oyn = length(unique(dataest1$GEOID))
oyInd = which(dataest5$GEOID[1:n] %in% dataest1$GEOID[1:oyn])

#standardize prescribing rate to use as covariate in W
#note prescribing rate only available for 2010 - 2019
pres.rate = rep(0,n*10)
for(i in 1:10){
  temp = data$PrescribingRate[((i+3-1)*n+1):((i+3)*n)]
  pres.rate[((i-1)*n+1):(i*n)] = (temp-mean(temp))/sd(temp)
}


### ACS 5-year estimates for 2009 - 2018
W5.obs = cbind(dataest5$Poverty, dataest5$Unemployed, dataest5$HS_Ed, dataest5$Food_Stamps)
W5.se = cbind(datamoe5$Poverty, datamoe5$Unemployed, datamoe5$HS_Ed, datamoe5$Food_Stamps)
W5.se = W5.se/1.645

### ACS 1-year estimates for 2009 - 2018
W1.obs = cbind(dataest1$Poverty, dataest1$Unemployed, dataest1$HS_Ed, dataest1$Food_Stamps)
W1.se = cbind(datamoe1$Poverty, datamoe1$Unemployed, datamoe1$HS_Ed, datamoe1$Food_Stamps)
W1.se = W1.se/1.645


### load in the values from the SAMSHA survey data
ell.lb=c(-3,1,5,9,10,11,12) ### first year of survey
ell.ub=c(0,4,8,10,11,12,13) ### last year of survey
ell.rate = c(-6/4,10/4,26/4,19/2,21/2,23/2,25/2)
S = c(.05,.055,0.052,.047,0.051,0.041,0.043) ### survey estimates
S.se=c(.0025,.0024,.0026,.0041,.0037,.0033,.0043) ### survey standard errors

y.T = data$Treatment ### standardized treatment counts due to DUA
y.D = data$Deaths
Pop = data$Population
Year = data$Year

#### identify the counties that have censored treatment values
#### using the simulsated data none are censored 
#### with disaggregated data censored by age group use commented out code
censC = rep(0,n*T)
censA = rep(0,n*T)
censBoth = rep(0,n*T)
censNone = rep(1,n*T)
cenlb = y.T
cenub = y.T
#censC = 1*(is.na(data$Treatment_u20))*(1-is.na(data$Treatment_o20))
#censA = 1*(1-is.na(data$Treatment_u20))*(is.na(data$Treatment_o20))
#censBoth = 1*(is.na(data$Treatment_u20))*(is.na(data$Treatment_o20))
#censNone = 1*(1-is.na(data$Treatment_u20))*(1-is.na(data$Treatment_o20))
#cenlb = rep(NA,length(y.T))
#cenub = rep(NA,length(y.T))
#cenlb[which(censC==1)]=data$Treatment_o20[which(censC==1)]
#cenub[which(censC==1)]=data$Treatment_o20[which(censC==1)]+9
#cenlb[which(censA==1)]=data$Treatment_u20[which(censA==1)]
#cenub[which(censA==1)]=data$Treatment_u20[which(censA==1)]+9
#cenlb[which(censBoth==1)]=1
#cenub[which(censBoth==1)]=18
#cenlb[which(is.na(cenlb))]=y.T[which(is.na(cenlb))]
#cenub[which(is.na(cenub))]=y.T[which(is.na(cenub))]

cen=censA+censC+censBoth
bd=cbind(cenlb,cenub)

Y = cbind(y.D,y.T)

#### load in the adjacency matrix for the Ohio counties
A<-read.csv(file = "OhioAdjacency.csv",header = FALSE)


### create design matrices for death and treatment rate among PWMO
X = matrix(0,n*T,T)
for(i in 1:T){
  X[((i-1)*n+1):((i-1)*n+n),i]=rep(1,n)
}
X.T = cbind(X,data$HPSA,data$MU)
X.D = cbind(X,1*(data$Interstates>0),data$HIDTA,data$MSA)


############################################################
###########           Set up            ############
############################################################

#### information needed for spatial structure in the CAR models
num<-colSums(A)
adj<-NULL
for(j in 1:n){
  adj<-c(adj,which(A[j,]==1))
}
adj<-as.vector(adj)
num<-as.vector(num)
weights<-1+0*num

############################################################
###########           WRite NIMBLE model       ############
############################################################

model_code=nimbleCode({
  for(t in 1:1){
    for(i in 1:R){
      pi.D[(t-1)*R+i] <- ilogit(inprod(X.D[(t-1)*R+i,],beta.D[])+(f[(t-1)*R+i]+mu.f[(t-1)*R+i])+eps.D[(t-1)*R+i])
      pi.T[(t-1)*R+i] <- ilogit(inprod(X.T[(t-1)*R+i,],beta.T[])+(g[(t-1)*R+i]+mu.g[(t-1)*R+i])+eps.T[(t-1)*R+i])
      mu.u[(t-1)*R+i] <- inprod(W.cent[(t+1)*R+i,1:(gp-1)],gamma[1:(gp-1)])
      lambda[(t-1)*R+i] <- exp((u[(t-1)*R+i]+mu.u[(t-1)*R+i])+v[(t-1)*R+i])
      y1[(t-1)*R+i] ~ dbinom(pi.D[(t-1)*R+i],N[(t-1)*R+i])
      cen[(t-1)*R+i] ~ dinterval(y2[(t-1)*R+i],bd[(t-1)*R+i,1:2])
      y2[(t-1)*R+i] ~ dbinom(pi.T[(t-1)*R+i],N[(t-1)*R+i])
      N[(t-1)*R+i] ~ dbinom(mu[t]*lambda[(t-1)*R+i],P[(t-1)*R+i])
      v[(t-1)*R+i] ~ dnorm(0,tau=tau.v)
      eps.D[(t-1)*R+i] ~ dnorm(0,tau=tau.D)
      eps.T[(t-1)*R+i] ~ dnorm(0,tau=tau.T)
      mu.f[(t-1)*R+i]<-0
      mu.g[(t-1)*R+i]<-0
    }
    mu[t] <- beta[1]+beta[2]*t
    f[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.f,zero_mean=1)
    u[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.u,zero_mean=1)
    g[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.g,zero_mean=1)
  }
  
  ### years 2 and 3 need to be different since OARRS data starts at 2010
  for(t in 2:3){
    for(i in 1:R){
      pi.D[(t-1)*R+i] <- ilogit(inprod(X.D[(t-1)*R+i,],beta.D[])+(f[(t-1)*R+i]+mu.f[(t-1)*R+i])+eps.D[(t-1)*R+i])
      pi.T[(t-1)*R+i] <- ilogit(inprod(X.T[(t-1)*R+i,],beta.T[])+(g[(t-1)*R+i]+mu.g[(t-1)*R+i])+eps.T[(t-1)*R+i])
      lambda[(t-1)*R+i] <- exp((u[(t-1)*R+i]+mu.u[(t-1)*R+i])+v[(t-1)*R+i])
      y1[(t-1)*R+i] ~ dbinom(pi.D[(t-1)*R+i],N[(t-1)*R+i])
      cen[(t-1)*R+i] ~ dinterval(y2[(t-1)*R+i],bd[(t-1)*R+i,1:2])
      y2[(t-1)*R+i] ~ dbinom(pi.T[(t-1)*R+i],N[(t-1)*R+i])
      N[(t-1)*R+i] ~ dbinom(mu[t]*lambda[(t-1)*R+i],P[(t-1)*R+i])
      v[(t-1)*R+i] ~ dnorm(0,tau=tau.v)
      eps.D[(t-1)*R+i] ~ dnorm(0,tau=tau.D)
      eps.T[(t-1)*R+i] ~ dnorm(0,tau=tau.T)
      mu.f[(t-1)*R+i] <- phi.f*f[(t-2)*R+i]
      mu.g[(t-1)*R+i] <- phi.g*g[(t-2)*R+i]
      mu.u[(t-1)*R+i] <- inprod(W.cent[(t+1)*R+i,1:(gp-1)],gamma[1:(gp-1)])+phi.u*(u[(t-2)*R+i]-inprod(W.cent[(t-0)*R+i,1:(gp-1)],gamma[1:(gp-1)]))
    }
    mu[t] <- beta[1]+beta[2]*t
    f[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.f,zero_mean=1)
    u[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.u,zero_mean=1)
    g[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.g,zero_mean=1)
  }
  
  ### year 4 slightly different in what we subtract for u
  for(t in 4:4){
    for(i in 1:R){
      pi.D[(t-1)*R+i] <- ilogit(inprod(X.D[(t-1)*R+i,],beta.D[])+(f[(t-1)*R+i]+mu.f[(t-1)*R+i])+eps.D[(t-1)*R+i])
      pi.T[(t-1)*R+i] <- ilogit(inprod(X.T[(t-1)*R+i,],beta.T[])+(g[(t-1)*R+i]+mu.g[(t-1)*R+i])+eps.T[(t-1)*R+i])
      lambda[(t-1)*R+i] <- exp((u[(t-1)*R+i]+mu.u[(t-1)*R+i])+v[(t-1)*R+i])
      y1[(t-1)*R+i] ~ dbinom(pi.D[(t-1)*R+i],N[(t-1)*R+i])
      cen[(t-1)*R+i] ~ dinterval(y2[(t-1)*R+i],bd[(t-1)*R+i,1:2])
      y2[(t-1)*R+i] ~ dbinom(pi.T[(t-1)*R+i],N[(t-1)*R+i])
      N[(t-1)*R+i] ~ dbinom(mu[t]*lambda[(t-1)*R+i],P[(t-1)*R+i])
      v[(t-1)*R+i] ~ dnorm(0,tau=tau.v)
      eps.D[(t-1)*R+i] ~ dnorm(0,tau=tau.D)
      eps.T[(t-1)*R+i] ~ dnorm(0,tau=tau.T)
      mu.f[(t-1)*R+i] <- phi.f*f[(t-2)*R+i]
      mu.g[(t-1)*R+i] <- phi.g*g[(t-2)*R+i]
      mu.u[(t-1)*R+i] <- inprod(W.cent[(t+1)*R+i,1:(gp-1)],gamma[1:(gp-1)])+pres.rate[(t-4)*R+i]*gamma[gp]+phi.u*(u[(t-2)*R+i]-inprod(W.cent[(t-0)*R+i,1:(gp-1)],gamma[1:(gp-1)]))
    }
    mu[t] <- beta[1]+beta[2]*t
    f[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.f,zero_mean=1)
    u[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.u,zero_mean=1)
    g[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.g,zero_mean=1)
  }
  
  
  for(t in 5:T){
    for(i in 1:R){
      pi.D[(t-1)*R+i] <- ilogit(inprod(X.D[(t-1)*R+i,],beta.D[])+(f[(t-1)*R+i]+mu.f[(t-1)*R+i])+eps.D[(t-1)*R+i])
      pi.T[(t-1)*R+i] <- ilogit(inprod(X.T[(t-1)*R+i,],beta.T[])+(g[(t-1)*R+i]+mu.g[(t-1)*R+i])+eps.T[(t-1)*R+i])
      lambda[(t-1)*R+i] <- exp((u[(t-1)*R+i]+mu.u[(t-1)*R+i])+v[(t-1)*R+i])
      y1[(t-1)*R+i] ~ dbinom(pi.D[(t-1)*R+i],N[(t-1)*R+i])
      cen[(t-1)*R+i] ~ dinterval(y2[(t-1)*R+i],bd[(t-1)*R+i,1:2])
      y2[(t-1)*R+i] ~ dbinom(pi.T[(t-1)*R+i],N[(t-1)*R+i])
      N[(t-1)*R+i] ~ dbinom(mu[t]*lambda[(t-1)*R+i],P[(t-1)*R+i])
      v[(t-1)*R+i] ~ dnorm(0,tau=tau.v)
      eps.D[(t-1)*R+i] ~ dnorm(0,tau=tau.D)
      eps.T[(t-1)*R+i] ~ dnorm(0,tau=tau.T)
      mu.f[(t-1)*R+i] <- phi.f*f[(t-2)*R+i]
      mu.g[(t-1)*R+i] <- phi.g*g[(t-2)*R+i]
      mu.u[(t-1)*R+i] <- inprod(W.cent[(t+1)*R+i,1:(gp-1)],gamma[1:(gp-1)])+pres.rate[(t-4)*R+i]*gamma[gp]+phi.u*(u[(t-2)*R+i]-inprod(W.cent[(t-0)*R+i,1:(gp-1)],gamma[1:(gp-1)])-pres.rate[(t-5)*R+i]*gamma[gp])
    }
    mu[t] <- beta[1]+beta[2]*t
    f[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.f,zero_mean=1)
    u[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.u,zero_mean=1)
    g[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.g,zero_mean=1)
  }
  
  #### need to add in model for W
  for(t in 3:T){
    for(j in 1:(gp-1)){ ### gp = number of variables in W
      ### need to make W go back to 2005 then cut off at 2007 above; Wobs starts at 2009
      for(i in 1:R){
        mu.Wobs[(t-3)*R+i,j] <- 1/5*(W[(t-3)*R+i,j]+W[(t-2)*R+i,j]+W[(t-1)*R+i,j]+W[(t-0)*R+i,j]+W[(t+1)*R+i,j])
        W5obs[(t-3)*R+i,j] ~ dnorm(mu.Wobs[(t-3)*R+i,j], sd=W5se[(t-3)*R+i,j])
      }
      for(l in 1:oyn){
        b[l,j,t] <- oyInd[l]
        W1obs[(t-3)*oyn+l,j] ~ dnorm(W[(t+1)*R+b[l,j,t],j],sd=W1se[(t-3)*oyn+l,j])
      }
    }
  }
  for(t in 1:(T+2)){
    for(j in 1:(gp-1)){
      W.cent[((t-1)*R+1):((t-1)*R+R),j] <- ( W[((t-1)*R+1):((t-1)*R+R),j]-mean( W[((t-1)*R+1):((t-1)*R+R),j]))/sd( W[((t-1)*R+1):((t-1)*R+R),j])
    }
  }
  for(j in 1:(gp-1)){
    for(t in 1:(T+2)){
      for(i in 1:R){
        W[(t-1)*R+i,j] ~ T(dnorm(mu.W[t,j],tau=tau.W[i,j]),0,100)
      }
      mu.W[t,j] ~ dunif(0,100)
    }
  }
  
  for(l in 1:L){
    S[l] ~ T(dnorm(beta[1]+beta[2]*ell.rate[l],sd=S.se[l]),0,1)
  }
  for(i in 1:bpT){
    beta.T[i] ~ dflat()
  }
  for(i in 1:bpD){
    beta.D[i] ~ dflat()
  }
  for(i in 1:2){
    beta[i] ~ dflat()
  }
  for(i in 1:(gp-1)){
    gamma[i] ~ dflat()
    for(j in 1:R){
      tau.W[j, i] ~ dgamma(.5, .5)
    }
  }
  gamma[gp] ~ dflat()
  tau.v ~ dgamma(.5,.5)
  tau.u ~ dgamma(.5,.5)
  tau.f ~ dgamma(.5,.5)
  tau.g ~ dgamma(.5,.5)
  tau.D ~ dgamma(.5,.5)
  tau.T ~ dgamma(.5,.5)
  phi.f ~ dunif(0,1)
  phi.u ~ dunif(0,1)
  phi.g ~ dunif(0,1)
})


mod_constants=list(bd=bd,gp = (length(W5.obs[1,])+1),bpT=length(X.T[1,]),bpD = length(X.D[1,]),R=n,T=T,L=length(ell.rate),ell.rate=ell.rate,P=Pop,adj=adj,X.T=X.T,X.D=X.D,num=num,W5se=W5.se,W1se=W1.se,S.se=S.se,oyn=oyn,oyInd=oyInd)
mod_data=list(y1=Y[,1],y2=Y[,2],cen=cen,S=S,W5obs=W5.obs,W1obs=W1.obs,pres.rate=pres.rate)

##### set the initial values

### need initial values for the censored treatment counts; this will be ignored if not censored
yTinit = y.T
yTinit[which(censC==1)] = data$Treatment_o20[which(censC==1)]+5
yTinit[which(censA==1)] = data$Treatment_u20[which(censA==1)]+5
yTinit[which(censBoth==1)] = 10


Ninit = floor(1/2*(5*yTinit+10*Y[,1])+100) ### check that N is not greater than P and also not less than T or D
II=which(Ninit>Pop)
if(length(II)>0){
  Ninit[II]=floor(.1*Pop[II])
}
II = which(Ninit<Y[,1])
if(length(II)>0){
  Ninit[II] = Y[II,1]+100
}
II = which(Ninit<Y[,2])
if(length(II)>0){
  Ninit[II] = Y[II,2]+100
}

ginit = rep(0,n*T)
finit = rep(0,n*T)
uinit = rep(0,n*T)


# Set initial values.
mod_inits=list(N=Ninit,beta.T=c(rep(-3.5,T),rep(0,length(X.T[1,])-T)),beta.D=c(rep(-6,T),rep(0,length(X.D[1,])-T)),beta=c(.05,-.001),tau.v=.1,tau.f=.5,tau.g=.5,tau.u=.5,
               u=uinit,v=rep(0,n*T),f=finit,g=ginit,tau.D=.1,tau.T=.1,eps.T=rep(0,n*T),eps.D=rep(0,n*T),W = rbind(W5.obs[1:(4*88),],W5.obs),mu.W = matrix(50,T+2,4),tau.W = matrix(1,n,length(W5.obs[1,])),
               gamma=rep(0,(length(W5.obs[1,])+1)),phi.u=.5,phi.f=.5,phi.g=.5,y2=floor(.5*(bd[,1]+bd[,2])))


# Build the model.
nimble_model <- nimbleModel(model_code, mod_constants,mod_data,mod_inits)
compiled_model <- compileNimble(nimble_model,resetFunctions = TRUE)

# Set up samplers.
mcmc_conf <- configureMCMC(nimble_model,monitors=c('beta.T','beta.D','beta','N','u','f','g','gamma','gamma.pr','W','mu.W','y2','pi.T','pi.D','lambda'),control=list(adaptive=TRUE,adaptInterval=100,sliceWidths=5,sliceMaxSteps=1000,maxContractions=10000000,sliceAdaptWidthMaxIter=0,sliceAdaptFactorMaxIter=0),useConjugacy = TRUE)

### change to AF Slice sampler for N_{i1},...,N_{iT}
for(i in 1:n){
  Nvec = rep(0,T)
  for(j in 1:T){
    Nvec[j]=paste("N[", toString(n*(j-1)+i), "]",sep="")
  }
  mcmc_conf$removeSamplers(Nvec)
  mcmc_conf$addSampler(target=Nvec,type='AF_slice',control=list(adaptive=TRUE,sliceWidths=rep(100,T),sliceMaxSteps=10000,maxContractions=1000000))
}


nimble_mcmc<-buildMCMC(mcmc_conf)
compiled_mcmc<-compileNimble(nimble_mcmc, project = nimble_model,resetFunctions = TRUE)


# Run the model 
MCS=1000000
st = Sys.time()
samples=runMCMC(compiled_mcmc,inits=mod_inits,
                   nchains = 1, nburnin=MCS/2,niter = MCS,samplesAsCodaMCMC = TRUE,thin=50,
                   summary = FALSE, WAIC = FALSE,progressBar=TRUE) 
Sys.time()-st


save(samples,file="OhioOutput.Rda")

