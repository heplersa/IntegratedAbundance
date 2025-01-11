# LOAD R PACKAGES
library(tidyverse) # data manipulation and visualization
library(nimble) # Bayesian inference for multi-level models

# IMPORT PRE-PROCESSED DATA USED TO FIT MODEL (data_for_analysis.Rda is generated from the WA_abundance_data.R file)
load("WAprevalence/data/data_for_analysis.Rda")

# DEFINE NIMBLE MODEL
model_code <- nimbleCode({
  
  for(t in 1:1){ # first year
    
    for(i in 1:R){ # counties
      
      for(j in 1:K){ # outcomes
        
        y[(t-1)*R+i,j] ~ dbinom(pi[(t-1)*R+i,j], N[(t-1)*R+i])
        mu.f[(t-1)*R+i,j] <- 0
        
      }
      
      # outcome error term
      eps[(t-1)*R+i, 1:3] ~ dmnorm(mean = mean.eps[1:3], cov = cov.eps[1:3, 1:3])
      
      # logistic link for binomial outcome(s)
      pi[(t-1)*R+i,1] <- ilogit(beta[t,1] + (f[(t-1)*R+i,1] + mu.f[(t-1)*R+i,1]) + eps[(t-1)*R+i,1])
      pi[(t-1)*R+i,2] <- ilogit(beta[t,2] + (f[(t-1)*R+i,2] + mu.f[(t-1)*R+i,2]) + eps[(t-1)*R+i,2])
      pi[(t-1)*R+i,3] <- ilogit(beta[t,3] + (f[(t-1)*R+i,3] + mu.f[(t-1)*R+i,3]) + eps[(t-1)*R+i,3])
      
      # latent counts (process model)
      mu.u[(t-1)*R+i] <- 0
      lambda[(t-1)*R+i] <- exp((u[(t-1)*R+i] + mu.u[(t-1)*R+i]) + v[(t-1)*R+i])
      N[(t-1)*R+i] ~ dbinom(mu[t]*lambda[(t-1)*R+i],P[(t-1)*R+i])
      v[(t-1)*R+i] ~ dnorm(0,tau.v)
      
      
    }
    
    # mean state-wide average risk of misuse in year t
    mu[t] <- ilogit(beta.mu[1] + beta.mu[2]*t)
    
  }
  
  for(t in 2:T){ # remaining years
    
    for(i in 1:R){ # counties
      
      for(j in 1:K){ # outcomes
        
        y[(t-1)*R+i,j] ~ dbinom(pi[(t-1)*R+i,j], N[(t-1)*R+i])
        mu.f[(t-1)*R+i,j] <- phi.f[j]*f[(t-2)*R+i, j]
        
      }
      
      # outcome error term
      eps[(t-1)*R+i, 1:3] ~ dmnorm(mean = mean.eps[1:3], cov = cov.eps[1:3, 1:3])
      
      # logistic link for binomial outcome(s)
      pi[(t-1)*R+i,1] <- ilogit(beta[t,1] + (f[(t-1)*R+i,1] + mu.f[(t-1)*R+i,1]) + eps[(t-1)*R+i,1])
      pi[(t-1)*R+i,2] <- ilogit(beta[t,2] + (f[(t-1)*R+i,2] + mu.f[(t-1)*R+i,2]) + eps[(t-1)*R+i,2])
      pi[(t-1)*R+i,3] <- ilogit(beta[t,3] + (f[(t-1)*R+i,3] + mu.f[(t-1)*R+i,3]) + eps[(t-1)*R+i,3])
      
      # latent counts (process model)
      mu.u[(t-1)*R+i] <- phi.u*u[(t-2)*R+i] 
      lambda[(t-1)*R+i] <- exp((u[(t-1)*R+i] + mu.u[(t-1)*R+i]) + v[(t-1)*R+i])
      N[(t-1)*R+i] ~ dbinom(mu[t]*lambda[(t-1)*R+i], P[(t-1)*R+i])
      v[(t-1)*R+i] ~ dnorm(0, tau.v)
      
    }
    
    # mean state-wide risk of misuse in year t
    mu[t] <- ilogit(beta.mu[1] + beta.mu[2]*t)
    
  }
  
  # state-wide survey data model
  for(l in 1:L){
    
    S[l] ~ dnorm(beta.mu[1]+beta.mu[2]*ell.rate[l], sd=S.se[l])
    
  }
  
  # spatial random effects for data level
  for(j in 1:K){
    
    for(t in 1:T){
      
      f[((t-1)*R+1):((t-1)*R+R), j] ~ dcar_normal(adj=adj[], 
                                                  num=num[],
                                                  tau=tau.f[j],
                                                  zero_mean=1)
      
    }
    
  }
  
  # spatial random effect for process level
  for(t in 1:T){
    
    u[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], 
                                             num=num[],
                                             tau=tau.u,
                                             zero_mean=1)
    
  }
  
  # time-varying intercepts for data level
  for(j in 1:K){
    
    for(t in 1:T){
      
      beta[t, j] ~ dflat() 
      
    }
    
    tau[j] ~ dgamma(.5,.5) # variance parameters for unexplained heterogeneity in data level
    phi.f[j] ~ dunif(0,1) # auto-regressive parameters for spatial random effect in data level
    tau.f[j] ~ dgamma(.5,.5) # variance parameters for for spatial random effect in data level
    
  }
  
  beta.mu[1:2] ~ dmnorm(mean = mean.mu[1:2], cov = cov.mu[1:2, 1:2])
  
  cov.eps[1:3, 1:3] ~ dwish(R=cov.eps.R[1:3, 1:3], df=3)
  
  tau.u ~ dgamma(.5,.5) # variance parameter for for spatial random effect in process level
  phi.u ~ dunif(0,1) # auto-regressive parameter for spatial random effect in process level
  tau.v ~ dgamma(.5,.5) # variance parameter for unexplained heterogeneity in process level
  
})

# DEFINE NIMBLE CONSTANTS, DATA, and INITS
n <- length(num) # number of WA counties
T <- length(2017:2022) # number of years
K <-  dim(yfit[, 3:5])[2] # number of outcomes

mod_constants <- list(R = n,
                      T = T,
                      K = K,
                      L = length(ell.rate),
                      ell.rate = ell.rate,
                      P = yfit$pop,
                      num = num,
                      adj = adj, 
                      S.se = logit_S.se,
                      mean.mu = rep(0, 2),
                      cov.mu = 10^4*diag(2),
                      mean.eps = rep(0, 3),
                      cov.eps.R = diag(3)
)

mod_data <- list(y=as.matrix(yfit[,3:5]),
                 S=logit_S
)

# specify initial values
logit_beta.mu.init <- lm(logit_S~ell.rate)$coefficients
beta.init <- matrix(data = 0, nrow = T, ncol = K)
Ninit <- floor(yfit$pop*ilogit(logit_beta.mu.init[1]))
finit <- matrix(data = 0, nrow = n*T, ncol = K)
uinit <- rep(0, n*T)

# constraints on domain of various outcomes. Make sure inits for N are less than P and greater than the outcome counts
II <- which(Ninit > yfit$pop)

if(length(II) > 0){
  
  Ninit[II] <- floor(.1*yfit[II, 6])
  
}

II <- which(Ninit < yfit[,3])

if(length(II) > 0){
  
  Ninit[II] <- yfit[II,1] + 100
  
}

II <- which(Ninit < yfit[,4])

if(length(II) > 0){
  
  Ninit[II] <- yfit[II,2] + 100
  
}

II <- which(Ninit < yfit[,5])

if(length(II) > 0){
  
  Ninit[II] <- yfit[II,3] + 100
  
}

# set initial values.
mod_inits <- list(N = Ninit,
                  beta = beta.init,
                  cov.eps = diag(3),
                  tau.v = .1, 
                  u = uinit,
                  f = finit,
                  v = rep(0, n*T),
                  tau = rep(.1, K), 
                  phi.u = .5, 
                  phi.f = rep(.5, K), 
                  beta.mu = logit_beta.mu.init
)

# BUILD AND RUN NIMBLE MODEL
nimble_model <- nimbleModel(model_code,
                            mod_constants,
                            mod_data,
                            mod_inits)

compiled_model <- compileNimble(nimble_model,
                                resetFunctions = TRUE)

# Set up samplers.
mcmc_conf <- configureMCMC(nimble_model,
                           monitors=c('tau',
                                      'tau.u',
                                      'tau.f',
                                      'eps',
                                      'cov.eps',
                                      'v',
                                      'beta',
                                      'beta.mu',
                                      'N',
                                      'u',
                                      'f',
                                      'pi',
                                      'lambda',
                                      'mu'
                           ),
                           control=list(adaptive=TRUE,
                                        adaptInterval=100,
                                        sliceWidths=5,
                                        sliceMaxSteps=1000,
                                        maxContractions=10000000,
                                        sliceAdaptWidthMaxIter=0,
                                        sliceAdaptFactorMaxIter=0),
                           useConjugacy = TRUE)

# change to AF Slice sampler for N_{i1},...,N_{iT}
for(i in 1:n){
  
  Nvec <- rep(0,T)
  
  for(j in 1:T){
    
    Nvec[j] <- paste("N[", toString(n*(j-1)+i), "]",sep="")
    
  }
  
  mcmc_conf$removeSamplers(Nvec)
  mcmc_conf$addSampler(target=Nvec,
                       type='AF_slice',
                       control=list(adaptive=TRUE,
                                    sliceWidths=rep(100,T),
                                    sliceMaxSteps=10000,
                                    maxContractions=1000000)
  )
}


nimble_mcmc <- buildMCMC(mcmc_conf)
compiled_mcmc <- compileNimble(nimble_mcmc, project = nimble_model, resetFunctions = TRUE)


# Run the model 
set.seed(2)
MCS <- 1*10^6
st  <- Sys.time()
samples <- runMCMC(compiled_mcmc,
                   inits = mod_inits,
                   nchains = 1, 
                   nburnin=MCS/2,
                   niter = MCS,
                   samplesAsCodaMCMC = TRUE,
                   thin = 50,
                   summary = FALSE, 
                   WAIC = FALSE,
                   progressBar = TRUE,
                   setSeed = 2) 

Sys.time()-st

MCMCvis::MCMCtrace(samples, pdf = F, params = "cov.eps")

save(samples, file = "WAprevalence/output/MCMC_no_covariates.Rda")
