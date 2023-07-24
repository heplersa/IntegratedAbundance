
###### Implement Integrated Abundance Model for New York State Counties
###### Outcome data begins in 2007; ACS 1-year estimates begin in 2006, ACS 5-year estimates in 2009

library(nimble) # For MCMC computation using NIMBLE.
library(coda) # For manipulation of MCMC results.

#### note rxrate and homeless rate have some missing values and are modeled

############################################################
###########             Load data               ############
############################################################

load('data2run.Rda')
##### This loads the following variables:
#### y : nT x K matrix of the K outcome variables for each county and year
#### X, X.M, X.B, X.D, X.E : design matrix of covariate effects that are *not* ACS estimates for misuse, methadone, buprenorphine, death, and ED visits
#### W5.obs, W5.se : matrix of ACS 5-year estimates and standard errors for each county from 2009 - 2018
#### W1.obs, W1.se : matrix of ACS 1-year estimates and standard errors for the counties that have these values from 2006 - 2018
#### W1.state, W1.state.se : matrix of ACS state-level 1-year estimates and standard errors for 2006 - 2018
#### A : adjacency matrix indicating what counties share a border
#### oyInd, oyInd.nhw, oyInd.vet : matrix that indicates which counties have 1-year ACS estimates for each year. Separate matrices are used for non-Hispanic White and percent veteran as they had additional missing values
#### ell.rate is the vector of the trend rate in the survey data
#### S, S.se are the survey estimates and standard errors
#### Pop is a vector of the population for each county in each year

wp = ncol(W5.obs)-1 ### number of ACS predictor variables in misuse since one ACS variable used in an outcome only

T = dim(y)[1]/dim(A)[1] ### data runs from 2007 - 2018; 1-year ACS data 2006-2018, 5-year ACS 2009-2018
#### need to model latent ACS variables from 2005 - 2018

### A is the adjacency matrix indicating neighboring counties
num<-colSums(A)
adj<-NULL
for(j in 1:n){
  adj<-c(adj,which(A[j,]==1))
}
adj<-as.vector(adj)
num<-as.vector(num)
weights<-1+0*num

### map what rows of the observed 1-year estimates in W1 go to what rows of the latent yearly W
#### note W starts in 2005 but W1.obs starts in 2006
oyn = rowSums(1-is.na(oyInd))
W1map = n+oyInd[1,1:oyn[1]] ### has what subset of rows from W correspond to data in W1
for(t in 2:(T+1)){
  W1map = c(W1map,(t)*n+oyInd[t,1:oyn[t]])
}
oyn.vet = rowSums(1-is.na(oyInd.vet))
W1map.vet = n+oyInd.vet[1,1:oyn.vet[1]] ### has what subset of rows from W with no missing data on pctveteran
for(t in 2:(T+1)){
  W1map.vet = c(W1map.vet,(t)*n+oyInd.vet[t,1:oyn.vet[t]])
}
### what rows of W1 are not missing for percent veteran variable
vetInd = which(is.na(W1.obs[,which(W.var=="pctveteran")])==FALSE)
oyn.nhw = rowSums(1-is.na(oyInd.nhw))
W1map.nhw = n+oyInd.nhw[1,1:oyn.nhw[1]] ### has what subset of rows from W with no missing data on pctnhweran
for(t in 2:(T+1)){
  W1map.nhw = c(W1map.nhw,(t)*n+oyInd.nhw[t,1:oyn.nhw[t]])
}
### what rows of W1 are not missing for percent non-Hispanic white variable
nhwInd = which(is.na(W1.obs[,which(W.var=="pctnhwhite")])==FALSE)


############################################################
###########           WRite NIMBLE model       ############
############################################################

model_code=nimbleCode({
 
  for(t in 1:1){
    for(i in 1:R){
      for(j in 1:(K-1)){
        y[(t-1)*R+i,j] ~ dbinom(pi[(t-1)*R+i,j],N[(t-1)*R+i])
        mu.f[(t-1)*R+i,j]<-0
        eps[(t-1)*R+i,j] ~ dnorm(0,tau[j])
      }
      y[(t-1)*R+i,K] ~ dpois(N[(t-1)*R+i]*pi[(t-1)*R+i,K])
      mu.f[(t-1)*R+i,K]<-0
      eps[(t-1)*R+i,K] ~ dnorm(0,tau[K])
      pi[(t-1)*R+i,1] <- ilogit(beta[t,1]+X.M[(t-1)*R+i]*beta.M+(f[(t-1)*R+i,1]+mu.f[(t-1)*R+i,1])+eps[(t-1)*R+i,1])
      pi[(t-1)*R+i,2] <- ilogit(beta[t,2]+X.B[(t-1)*R+i]*beta.B+(f[(t-1)*R+i,2]+mu.f[(t-1)*R+i,2])+eps[(t-1)*R+i,2])
      pi[(t-1)*R+i,3] <- ilogit(beta[t,3]+inprod(X.D[(t-1)*R+i,1:dp],beta.D[1:dp])+W.std[(t+1)*R+i,(wp+1)]*alpha.D+(f[(t-1)*R+i,3]+mu.f[(t-1)*R+i,3])+eps[(t-1)*R+i,3])
      log(pi[(t-1)*R+i,K]) <- beta[t,K]+inprod(X.E[(t-1)*R+i,1:ep],beta.E[1:ep])+W.std[(t+1)*R+i,(wp+1)]*alpha.E+(f[(t-1)*R+i,K]+mu.f[(t-1)*R+i,K])+eps[(t-1)*R+i,K]
      mu.u[(t-1)*R+i] <- inprod(X[(t-1)*R+i,1:gp],beta.U[1:gp])+inprod(W.std[(t+1)*R+i,1:wp],alpha.U[1:wp])
      lambda[(t-1)*R+i] <- exp((u[(t-1)*R+i]+mu.u[(t-1)*R+i])+v[(t-1)*R+i])
      N[(t-1)*R+i] ~ dbinom(mu[t]*lambda[(t-1)*R+i],P[(t-1)*R+i])
      v[(t-1)*R+i] ~ dnorm(0,tau.v)
      #### model for X[,1] (rxrate)
      X[(t-1)*R+i,1] ~ dnorm(rx[t],sd=1)
      #### model for X[,2] (homelessrate)
      X[(t-1)*R+i,2] ~ dnorm(hr[t],sd=1)
    }
    mu[t] <- beta.mu[1]+beta.mu[2]*t
    rx[t] ~ dflat()
    hr[t] ~ dflat()
  }
  for(t in 2:T){
    for(i in 1:R){
      for(j in 1:(K-1)){
        y[(t-1)*R+i,j] ~ dbinom(pi[(t-1)*R+i,j],N[(t-1)*R+i])
        mu.f[(t-1)*R+i,j]<-phi.f[j]*f[(t-2)*R+i,j]
        eps[(t-1)*R+i,j] ~ dnorm(0,tau[j])
      }
      y[(t-1)*R+i,K] ~ dpois(N[(t-1)*R+i]*pi[(t-1)*R+i,K])
      mu.f[(t-1)*R+i,K]<-phi.f[K]*f[(t-2)*R+i,K]
      eps[(t-1)*R+i,K] ~ dnorm(0,tau[K])
      pi[(t-1)*R+i,1] <- ilogit(beta[t,1]+X.M[(t-1)*R+i]*beta.M+(f[(t-1)*R+i,1]+mu.f[(t-1)*R+i,1])+eps[(t-1)*R+i,1])
      pi[(t-1)*R+i,2] <- ilogit(beta[t,2]+X.B[(t-1)*R+i]*beta.B+(f[(t-1)*R+i,2]+mu.f[(t-1)*R+i,2])+eps[(t-1)*R+i,2])
      pi[(t-1)*R+i,3] <- ilogit(beta[t,3]+inprod(X.D[(t-1)*R+i,1:dp],beta.D[1:dp])+W.std[(t+1)*R+i,wp+1]*alpha.D+(f[(t-1)*R+i,3]+mu.f[(t-1)*R+i,3])+eps[(t-1)*R+i,3])
      log(pi[(t-1)*R+i,K]) <- beta[t,K]+inprod(X.E[(t-1)*R+i,1:ep],beta.E[1:ep])+W.std[(t+1)*R+i,wp+1]*alpha.E+(f[(t-1)*R+i,K]+mu.f[(t-1)*R+i,K])+eps[(t-1)*R+i,K]
      mu.u[(t-1)*R+i] <- inprod(X[(t-1)*R+i,1:gp],beta.U[1:gp])+inprod(W.std[(t+1)*R+i,1:wp],alpha.U[1:wp])+phi.u*(u[(t-2)*R+i]-inprod(X[(t-2)*R+i,1:gp],beta.U[1:gp])-inprod(W.std[(t-0)*R+i,1:wp],alpha.U[1:wp]))
      lambda[(t-1)*R+i] <- exp((u[(t-1)*R+i]+mu.u[(t-1)*R+i])+v[(t-1)*R+i])
      N[(t-1)*R+i] ~ dbinom(mu[t]*lambda[(t-1)*R+i],P[(t-1)*R+i])
      v[(t-1)*R+i] ~ dnorm(0,tau.v)
      #### model for X[,1] (rxrate)
      X[(t-1)*R+i,1] ~ dnorm(rx[t],sd=1)
      #### model for X[,2] (homelessrate)
      X[(t-1)*R+i,2] ~ dnorm(hr[t],sd=1)
    }
    mu[t] <- beta.mu[1]+beta.mu[2]*t
    rx[t] ~ dflat()
    hr[t] ~ dflat()
  }
  
  for(t in 1:(T+2)){
    for(j in 1:(wp+1)){
      W.std[((t-1)*R+1):((t-1)*R+R),j] <- ( W[((t-1)*R+1):((t-1)*R+R),j]-mean( W[((t-1)*R+1):((t-1)*R+R),j]))/sd( W[((t-1)*R+1):((t-1)*R+R),j])
    }
  }
  
  
  #### model for 5-year ACS data
  for(t in 3:T){
    for(j in 1:(wp+1)){ ### = number of variables in W
      ### need to make W go back to 2005; W5obs starts at 2009
      for(i in 1:R){
        mu.Wobs[(t-3)*R+i,j] <- 1/5*(W[(t-3)*R+i,j]+W[(t-2)*R+i,j]+W[(t-1)*R+i,j]+W[(t-0)*R+i,j]+W[(t+1)*R+i,j])
        W5obs[(t-3)*R+i,j] ~ T(dnorm(mu.Wobs[(t-3)*R+i,j], sd=W5se[(t-3)*R+i,j]),0,100)
      }
    }
  }

  #### model for 1-year ACS data excluding pctveteran and nhwhite
  for(j in 1:(5)){
    for(i in 1:n0){
      W1obs[i,j] ~ T(dnorm(W[W1map[i],j],sd=W1se[i,j] ),0,100)
    }
  }
  
  #### model for 1-year ACS data for pctveteran
  for(j in 6:6){
    for(i in 1:n0.vet){
      W1obs[vetInd[i],j] ~ T(dnorm(W[W1map.vet[i],j],sd=W1se[vetInd[i],j] ),0,100)
    }
  }
  
  #### model for 1-year ACS data for nhw
  for(j in 7:7){
    for(i in 1:n0.nhw){
      W1obs[nhwInd[i],j] ~ T(dnorm(W[W1map.nhw[i],j],sd=W1se[nhwInd[i],j] ),0,100)
    }
  }
  
  #### model for 1-year ACS data for median age
  for(j in 8:8){
    for(i in 1:n0){
      W1obs[i,j] ~ T(dnorm(W[W1map[i],j],sd=W1se[i,j] ),0,100)
    }
  }
 

  #### model for state level 1-year ACS data for socioeconomic variables
  for(t in 1:(T+1)){
    for(j in 1:2){
      W1.state[t,j] ~ T(dnorm(mu.W.s[t+1,j],sd=W1.state.se[t,j]),0,100) ### note mu.W.s defined from 2005-
    }
  }

  
  #### model for latent W
  #### variables 1&2 are socioeconomic...3:6 are demographic are relatively stable over time
  for(j in 1:2){ ### socioeconomic variables
    for(t in 1:(T+2)){
      for(i in 1:R){
        W[(t-1)*R+i,j] ~ T(dnorm(mu.W.s[t,j]+g[(t-1)*R+i,j],tau.W[j]),0,100)
      }
      mu.W.s[t,j] ~ T(dnorm(mu.W.s0[j],tau.W0[j]),0,100)
      #### make g ICAR for each t,j
      g[((t-1)*R+1):((t-1)*R+R),j]~dcar_normal(adj=adj[], num=num[], tau=tau.g[j],zero_mean=1)
    }
    tau.g[j] ~ dgamma(.5,.5)
    tau.W[j] ~ dgamma(.5,.5)
    tau.W0[j] ~ dgamma(.5,.5)
    mu.W.s0[j] ~ dunif(0,100)
  }
  for(j in 3:(wp+1)){ ### demographic variables 
    for(i in 1:R){
      for(t in 1:(T+2)){
        W[(t-1)*R+i,j] ~ T(dnorm(mu.W.d[i,j-2],tau.W[(j)]),0,100)
      }
      mu.W.d[i,j-2] ~ dunif(0,100)
    }
    tau.W[j] ~ dgamma(.5,.5)
  }

  for(t in 1:T){
    for(j in 1:K){
      f[((t-1)*R+1):((t-1)*R+R),j]~dcar_normal(adj=adj[], num=num[], tau=tau.f[j],zero_mean=1)
    }
  }

  for(t in 1:T){
    u[((t-1)*R+1):((t-1)*R+R)]~dcar_normal(adj=adj[], num=num[], tau=tau.u,zero_mean=1)
  }
  
  for(l in 1:L){
    S[l] ~ T(dnorm(beta.mu[1]+beta.mu[2]*ell.rate[l],sd=S.se[l]),0,1)
  }
  for(j in 1:K){
    for(t in 1:T){
      beta[t, j] ~ dflat()
    }
    tau[j] ~ dgamma(.5,.5)
    phi.f[j] ~ dunif(0,1)
    tau.f[j] ~ dgamma(.5,.5)
  }
  
  beta.M ~ dflat()
  beta.B ~ dflat()
  
  for(p in 1:dp){
    beta.D[p] ~ dflat()
  }
  for(p in 1:ep){
    beta.E[p] ~ dflat()
  }
  for(p in 1:gp){
    beta.U[p] ~ dflat()
  }
  alpha.D ~ dflat()
  alpha.E ~ dflat()
  for(p in 1:wp){
    alpha.U[p] ~ dflat()
  }


  beta.mu[1] ~ dflat()
  beta.mu[2] ~ dflat()
  tau.u ~ dgamma(.5,.5)
  phi.u ~ dunif(0,1)
  tau.v ~ dgamma(.5,.5)
})


##### set the initial values
K = 4

y.temp = y
y.temp[is.na(y.temp[,2]),2]=0
beta.mu.init = lm(S~ell.rate)$coefficients
Ninit = floor(Pop*beta.mu.init[1])
 
II=which(Ninit>(.2*Pop))
if(length(II)>0){
  Ninit[II]=floor(.1*Pop[II])
}
II = which(Ninit<y[,1])
if(length(II)>0){
  Ninit[II] = y[II,1]+100
}
II = which(Ninit<y.temp[,2])
if(length(II)>0){
  Ninit[II] = y.temp[II,2]+100
}

### initial values for X and W
Xinit = X
Xinit[which(is.na(Xinit[,1])),1]=mean(X[,1],na.rm=TRUE)
Xinit[which(is.na(Xinit[,2])),2]=mean(X[,2],na.rm=TRUE)

###
Winit = matrix(0,n*(T+2),ncol(W5.obs))
Winit[1:(n*4),]=W5.obs[1:(n*4),]
Winit[(4*n+1):(n*(T+2)),]=W5.obs



finit = matrix(0,n*T,K)
uinit = rep(0,n*T)
beta.init = matrix(-4,T,K)

mod_constants=list(wp = wp, Wlb=rep(0,ncol(W5.obs)),Wub = rep(100,ncol(W5.obs)),vetInd=vetInd,n0.vet = sum(oyn.vet),nhwInd=nhwInd,n0.nhw=sum(oyn.nhw), n0 = nrow(W1.obs),W1map.vet=W1map.vet,W1map.nhw=W1map.nhw,W1map=W1map,oynsum = c(0,cumsum(oyn)),WisPsi=diag(K),gp = ncol(X),bp = ncol(X.B),mp = ncol(X.M),ep = ncol(X.E),dp = ncol(X.D),R=n,T=T,L=length(ell.rate),ell.rate=ell.rate,P=Pop,cov.mat=diag(n),num=num,adj=adj, zeros = rep(0, n),S.se=S.se,K=K,oyn=oyn,oyInd=oyInd,W5se=W5.se,W1.state.se=W1.state.se)
mod_data=list(y=y,S=S,X=X,X.D=X.D,X.E=X.E,X.M=as.vector(X.M),X.B=as.vector(X.B),W5obs=W5.obs,W1obs=W1.obs,W1.state=W1.state,W1se=W1.se)


# Set initial values.

mod_inits=list(rx = rep(.5,T), hr = rep(15,T),X=Xinit,W=Winit, N=Ninit,beta = beta.init,tau.v=.1,beta.D=rep(0,ncol(X.D)),beta.E=rep(0,ncol(X.E)),beta.M=rep(0,ncol(X.M)),beta.B=rep(0,ncol(X.B)),
               u=uinit,v=rep(0,n*T),f=finit,tau=rep(.1,K),g = matrix(0,n*(T+2),2),
               phi.u=.5,phi.f=rep(.5,K),beta.mu=beta.mu.init,beta.U = rep(0,ncol(X)),
               alpha.D=rep(0,1),alpha.E=rep(0,1),alpha.U = rep(0,(ncol(W5.obs)-1)))




# Build the model.
nimble_model <- nimbleModel(model_code, mod_constants,mod_data,mod_inits)
compiled_model <- compileNimble(nimble_model,resetFunctions = TRUE)

# Set up samplers.
mcmc_conf <- configureMCMC(nimble_model,monitors=c('tau.u','tau.f','beta','beta.mu','beta.D','beta.M','beta.B','beta.E','beta.U','N','pi','lambda','alpha.E','alpha.D','alpha.U','W','mu.W.s','mu.W.d'),control=list(adaptive=TRUE,adaptInterval=100,sliceWidths=5,sliceMaxSteps=1000,maxContractions=10000000,sliceAdaptWidthMaxIter=0,sliceAdaptFactorMaxIter=0),useConjugacy = TRUE)

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
MCS=500000
MCB=250000
MCT = 50
st = Sys.time()
samples=runMCMC(compiled_mcmc,inits=mod_inits,
                nchains = 1, nburnin=MCB,niter = MCS,samplesAsCodaMCMC = TRUE,thin=MCT,
                summary = FALSE, WAIC = FALSE,progressBar=TRUE) 
Sys.time()-st


 save(samples,file="NYOutput.Rda")

q()

