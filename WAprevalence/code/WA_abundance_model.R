# FIT INTEGRATED ABUNDANCE MODEL USING PROCESSED WA DATA #
# BRIAN N. WHITE #
# 2025-01-18 #

# LOAD R PACKAGES
library(tidyverse) # data manipulation and visualization
library(nimble) # Bayesian inference for multi-level models

# IMPORT PRE-PROCESSED DATA USED TO FIT MODEL (data_for_analysis.Rda is generated from the WA_abundance_data.R file)
load("WAprevalence/data/data_for_analysis.Rda")

# DEFINE NIMBLE MODEL
model_code <- nimbleCode({
  
  for(t in 1:1){ # first year
    
    for(i in 1:R){ # counties
      
      for(j in 1:2){ # outcomes (each patient contributes only one count)
        
        y[(t-1)*R+i,j] ~ dbinom(pi[(t-1)*R+i,j], N[(t-1)*R+i])
        mu.f[(t-1)*R+i,j] <- 0
        
      }
      
      # outcomes (each patient can contribute more than one count)
      
        # shift baseline time to be 2019 for ED visit outcome as no data available for 2017-2018
        y[(t-1+2)*R+i,3] ~ dpois(pi[(t-1+2)*R+i,3]*N[(t-1+2)*R+i])
        mu.f[(t-1+2)*R+i,3] <- 0
          
        # hospitalization outcome has data for all years so no need to shift time-index
        y[(t-1)*R+i,4] ~ dpois(pi[(t-1)*R+i,4]*N[(t-1)*R+i])
        mu.f[(t-1)*R+i,4] <- 0
        
        # poisson outcomes are censored in [1,9]
        censored_ed[(t-1+2)*R+i] ~ dinterval(y[(t-1+2)*R+i,3], c_ed[(t-1+2)*R+i, 1:2])
        censored_hosp[(t-1)*R+i] ~ dinterval(y[(t-1)*R+i,4], c_hosp[(t-1)*R+i, 1:2])
      
      # outcome error term
      eps[(t-1)*R+i, 1:K] ~ dmnorm(mean = mean.eps[1:K], cov = cov.eps[1:K, 1:K])
      
      # logistic link for binomial outcome(s)
      pi[(t-1)*R+i,1] <- ilogit(beta[t,1] + (f[(t-1)*R+i,1] + mu.f[(t-1)*R+i,1]) + eps[(t-1)*R+i,1])
      pi[(t-1)*R+i,2] <- ilogit(beta[t,2] + (f[(t-1)*R+i,2] + mu.f[(t-1)*R+i,2]) + eps[(t-1)*R+i,2])
      
      # log link for poisson outcome(s)
      pi[(t-1+2)*R+i,3] <- exp(beta[t+2,3] + (f[(t-1+2)*R+i,3] + mu.f[(t-1+2)*R+i,3]) + eps[(t-1+2)*R+i,3])
      pi[(t-1)*R+i,4] <- exp(beta[t,4] + (f[(t-1)*R+i,4] + mu.f[(t-1)*R+i,4]) + eps[(t-1)*R+i,4])
      
      # latent counts (process model)
      mu.u[(t-1)*R+i] <- 0
      lambda[(t-1)*R+i] <- exp((u[(t-1)*R+i] + mu.u[(t-1)*R+i]) + v[(t-1)*R+i])
      N[(t-1)*R+i] ~ dbinom(mu[t]*lambda[(t-1)*R+i],P[(t-1)*R+i])
      v[(t-1)*R+i] ~ dnorm(0,tau.v)
      
      
    }
    
    # mean state-wide average risk of misuse in year t
    mu[t] <- ilogit(beta.mu[1] + beta.mu[2]*t)
    
  }
  
  for(t in 2:T){ # remaining years (for outcomes with data available over whole study period)
    
    for(i in 1:R){ # counties
      
      for(j in 1:2){ # outcomes (each patient contributes only one count)
        
        y[(t-1)*R+i,j] ~ dbinom(pi[(t-1)*R+i,j], N[(t-1)*R+i])
        mu.f[(t-1)*R+i,j] <- phi.f[j]*f[(t-2)*R+i, j]
        
      }
      
      # outcomes (each patient can contribute more than one count)
        
        # hospitalization outcome has data for all years so no need to shift time-index
        y[(t-1)*R+i,4] ~ dpois(pi[(t-1)*R+i,4]*N[(t-1)*R+i])
        mu.f[(t-1)*R+i,4] <- phi.f[4]*f[(t-2)*R+i, 4]
        
        # poisson outcomes are censored in [1,9]
        censored_hosp[(t-1)*R+i] ~ dinterval(y[(t-1)*R+i,4], c_hosp[(t-1)*R+i, 1:2])
      
      # outcome error term
      eps[(t-1)*R+i, 1:K] ~ dmnorm(mean = mean.eps[1:K], cov = cov.eps[1:K, 1:K])
      
      # logistic link for binomial outcome(s)
      pi[(t-1)*R+i,1] <- ilogit(beta[t,1] + (f[(t-1)*R+i,1] + mu.f[(t-1)*R+i,1]) + eps[(t-1)*R+i,1])
      pi[(t-1)*R+i,2] <- ilogit(beta[t,2] + (f[(t-1)*R+i,2] + mu.f[(t-1)*R+i,2]) + eps[(t-1)*R+i,2])
      
      # log link for poisson outcome(s)
      pi[(t-1)*R+i,4] <- exp(beta[t,4] + (f[(t-1)*R+i,4] + mu.f[(t-1)*R+i,4]) + eps[(t-1)*R+i,4])

      # latent counts (process model)
      mu.u[(t-1)*R+i] <- phi.u*u[(t-2)*R+i] 
      lambda[(t-1)*R+i] <- exp((u[(t-1)*R+i] + mu.u[(t-1)*R+i]) + v[(t-1)*R+i])
      N[(t-1)*R+i] ~ dbinom(mu[t]*lambda[(t-1)*R+i], P[(t-1)*R+i])
      v[(t-1)*R+i] ~ dnorm(0, tau.v)
      
    }
    
    # mean state-wide risk of misuse in year t
    mu[t] <- ilogit(beta.mu[1] + beta.mu[2]*t)
    
  }
  
  for(t in 2:(T-2)){ # remaining years (for outcomes with data missing for 2018-2019)
    
    for(i in 1:R){ # counties
      
      # outcomes (each patient can contribute more than one count)
      
      # shift baseline time to be 2019 for ED visit outcome as no data available for 2017-2018
      y[(t-1+2)*R+i,3] ~ dpois(pi[(t-1+2)*R+i,3]*N[(t-1+2)*R+i])
      mu.f[(t-1+2)*R+i, 3] <- phi.f[3]*f[(t-2+2)*R+i, 3]
      
      # poisson outcomes are censored in [1,9]
      censored_ed[(t-1+2)*R+i] ~ dinterval(y[(t-1+2)*R+i,3], c_ed[(t-1+2)*R+i, 1:2])
      
      # log link for poisson outcome(s)
      pi[(t-1+2)*R+i,3] <- exp(beta[t+2,3] + (f[(t-1+2)*R+i,3] + mu.f[(t-1+2)*R+i,3]) + eps[(t-1+2)*R+i,3])
      
    }
    
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
    
    phi.f[j] ~ dunif(0,1) # auto-regressive parameters for spatial random effect in data level
    tau.f[j] ~ dgamma(.5,.5) # variance parameters for for spatial random effect in data level
    
  }
  
  beta.mu[1:2] ~ dmnorm(mean = mean.mu[1:2], cov = cov.mu[1:2, 1:2])
  
  cov.eps[1:K, 1:K] ~ dwish(R=cov.eps.R[1:K, 1:K], df=K)
  
  tau.u ~ dgamma(.5,.5) # variance parameter for for spatial random effect in process level
  phi.u ~ dunif(0,1) # auto-regressive parameter for spatial random effect in process level
  tau.v ~ dgamma(.5,.5) # variance parameter for unexplained heterogeneity in process level
  
})

# DEFINE NIMBLE CONSTANTS, DATA, and INITS
n <- length(num) # number of WA counties
T <- length(2017:2023) # number of years
K <-  dim(yfit[, c("pmp", "death", "ed", "hosp")])[2] # number of outcomes

# prepare constraint objects
cens_index_ed <- which(is.na(yfit$ed))[-c(1:78)] # note that 2017, 2018 NA but not censored
cens_index_hosp <- which(is.na(yfit$hosp))
censored_ed <- rep(1, n*T)
censored_hosp <- rep(1, n*T)
yinit <- yfit[, c("pmp", "death", "ed", "hosp")]
yinit[cens_index_ed, 3] <-  2
yinit[cens_index_hosp, 4] <-  2
c_ed <- matrix(0, nrow=n*T, ncol=2)
c_hosp <- matrix(0, nrow=n*T, ncol=2)
c_ed[cens_index_ed, 1] <- 0
c_ed[cens_index_ed, 2] <- 9
c_ed[-cens_index_ed, 1] <- -Inf
c_ed[-cens_index_ed, 2] <- Inf
c_hosp[cens_index_hosp, 1] <- 0
c_hosp[cens_index_hosp, 2] <- 9
c_hosp[-cens_index_hosp, 1] <- -Inf
c_hosp[-cens_index_hosp, 2] <- Inf

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
                      mean.eps = rep(0, K),
                      cov.eps.R = diag(K),
                      c_ed = c_ed,
                      c_hosp = c_hosp
)

mod_data <- list(y=as.matrix(yfit[,c("pmp", "death", "ed", "hosp")]),
                 censored_ed = censored_ed,
                 censored_hosp = censored_hosp,
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
  
  Ninit[II] <- floor(.1*yfit[II, "pop"])
  
}

II <- which(Ninit < yfit[,"pmp"])

if(length(II) > 0){
  
  Ninit[II] <- yfit[II,"pmp"] + 100
  
}

II <- which(Ninit < yfit[,"death"])

if(length(II) > 0){
  
  Ninit[II] <- yfit[II,"death"] + 100
  
}

II <- which(Ninit < yfit[,"ed"])

if(length(II) > 0){
  
  Ninit[II] <- yfit[II,"ed"] + 100
  
}

II <- which(Ninit < yfit[,"hosp"])

if(length(II) > 0){
  
  Ninit[II] <- yfit[II,"hosp"] + 100
  
}

# set initial values.
mod_inits <- list(y = as.matrix(yinit),
                  N = Ninit,
                  beta = beta.init,
                  cov.eps = diag(K),
                  tau.v = .1, 
                  u = uinit,
                  f = finit,
                  v = rep(0, n*T),
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

# set up samplers.
mcmc_conf <- configureMCMC(nimble_model,
                           monitors=c('tau.u',
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
set.seed(2025)
MCS <- 1*10^5
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
save(samples, file = "WAprevalence/output/mcmc/MCMC_no_covariates_2025_12_13.Rda")