# LOAD R PACKAGES
library(tidyverse) # data manipulation and visualization
library(nimble) # Bayesian inference for multi-level models

# IMPORT PRE-PROCESSED DATA USED TO FIT MODEL (data_for_analysis.Rda is generated from the NC_abundance_data.R file)
load("../data/data_for_analysis.Rda")

# DEFINE NIMBLE MODEL
model_code <- nimbleCode({
  
  for(t in 1:1){ # first year
    
    for(i in 1:R){ # counties
      
      for(j in c(1, K)){ 
        
        # uncensored outcome(s)
        y[(t-1)*R+i,j] ~ dbinom(pi[(t-1)*R+i,j], N[(t-1)*R+i])
        mu.f[(t-1)*R+i,j] <- 0
        
      }
      
        # censored outcome (treatment program use outcome censors non-zero counts less than 5)
        censored[(t-1)*R+i] ~ dinterval(y[(t-1)*R+i,2], c[(t-1)*R+i, 1:2])
        y[(t-1)*R+i,2] ~ dbinom(pi[(t-1)*R+i,2], N[(t-1)*R+i])
        mu.f[(t-1)*R+i,2] <- 0

      # outcome error term
      eps[(t-1)*R+i, 1:3] ~ dmnorm(mean = mean.eps[1:3], cov = cov.eps[1:3, 1:3])
      
      # logistic link for binomial outcome(s)
      pi[(t-1)*R+i,1] <- ilogit(beta[t,1] + inprod(X.D[(t-1)*R+i, 1:dp], beta.D[1:dp]) + (f[(t-1)*R+i,1] + mu.f[(t-1)*R+i,1]) + eps[(t-1)*R+i,1])
      pi[(t-1)*R+i,2] <- ilogit(beta[t,2] + inprod(X.T[(t-1)*R+i, 1:tp], beta.T[1:tp]) + (f[(t-1)*R+i,2] + mu.f[(t-1)*R+i,2]) + eps[(t-1)*R+i,2])
      pi[(t-1)*R+i,3] <- ilogit(beta[t,3] + inprod(X.B[(t-1)*R+i, 1:bp], beta.B[1:bp]) + (f[(t-1)*R+i,3] + mu.f[(t-1)*R+i,3]) + eps[(t-1)*R+i,3])
      
      # latent counts (process model)
      mu.u[(t-1)*R+i] <- inprod(W.cent[(t+3)*R+i,1:gp], gamma[1:gp]) # t+3=2016 for W.cent when t=1
      lambda[(t-1)*R+i] <- exp((u[(t-1)*R+i] + mu.u[(t-1)*R+i]) + v[(t-1)*R+i])
      N[(t-1)*R+i] ~ dbinom(mu[t]*lambda[(t-1)*R+i],P[(t-1)*R+i])
      v[(t-1)*R+i] ~ dnorm(0,tau.v)
      
      
    }
    
    # mean state-wide average risk of misuse in year t
    mu[t] <- ilogit(beta.mu[1] + beta.mu[2]*t)
    
  }
  
  for(t in 2:T){ # remaining years
    
    for(i in 1:R){ # counties
      
      for(j in c(1, K)){ 
        
        # uncensored outcome(s)
        y[(t-1)*R+i,j] ~ dbinom(pi[(t-1)*R+i,j], N[(t-1)*R+i])
        mu.f[(t-1)*R+i,j] <- phi.f[j]*f[(t-2)*R+i, j]
        
      }

        # censored outcome (treatment program use outcome censors non-zero counts less than 5)
        censored[(t-1)*R+i] ~ dinterval(y[(t-1)*R+i,2], c[(t-1)*R+i, 1:2])
        y[(t-1)*R+i,2] ~ dbinom(pi[(t-1)*R+i,2], N[(t-1)*R+i])
        mu.f[(t-1)*R+i,2] <- phi.f[2]*f[(t-2)*R+i, 2]
      
      # outcome error term
      eps[(t-1)*R+i, 1:3] ~ dmnorm(mean = mean.eps[1:3], cov = cov.eps[1:3, 1:3])
      
      # logistic link for binomial outcome(s)
      pi[(t-1)*R+i,1] <- ilogit(beta[t,1] + inprod(X.D[(t-1)*R+i, 1:dp], beta.D[1:dp]) + (f[(t-1)*R+i,1] + mu.f[(t-1)*R+i,1]) + eps[(t-1)*R+i,1])
      pi[(t-1)*R+i,2] <- ilogit(beta[t,2] + inprod(X.T[(t-1)*R+i, 1:tp], beta.T[1:tp]) + (f[(t-1)*R+i,2] + mu.f[(t-1)*R+i,2]) + eps[(t-1)*R+i,2])
      pi[(t-1)*R+i,3] <- ilogit(beta[t,3] + inprod(X.B[(t-1)*R+i, 1:bp], beta.B[1:bp]) + (f[(t-1)*R+i,3] + mu.f[(t-1)*R+i,3]) + eps[(t-1)*R+i,3])
      
      # latent counts (process model)
      mu.u[(t-1)*R+i] <- inprod(W.cent[(t+3)*R+i,1:gp],gamma[1:gp]) + phi.u*u[(t-2)*R+i] # t+3==2016 for W.cent when t=1
      lambda[(t-1)*R+i] <- exp((u[(t-1)*R+i] + mu.u[(t-1)*R+i]) + v[(t-1)*R+i])
      N[(t-1)*R+i] ~ dbinom(mu[t]*lambda[(t-1)*R+i], P[(t-1)*R+i])
      v[(t-1)*R+i] ~ dnorm(0, tau.v)
      
    }
    
    # mean state-wide risk of misuse in year t
    mu[t] <- ilogit(beta.mu[1] + beta.mu[2]*t)
    
  }
  
  
  # standardize W
  for(t in 1:(T+4)){ # T+4=2021 for W
    
    for(j in 1:gp){
      
      W.cent[((t-1)*R+1):((t-1)*R+R), j] <- ( W_ilogit[((t-1)*R+1):((t-1)*R+R),j]-mean(W_ilogit[((t-1)*R+1):((t-1)*R+R),j]))/sd( W_ilogit[((t-1)*R+1):((t-1)*R+R),j])
      
    }
    
  }
  
  # model for 5-year county data
  for(t in 1:T){
    
    for(j in 1:gp){
      
      # gp = number of variables in W
      # W, W1obs start in 2012; W5obs starts in 2016
      
      for(i in 1:R){
        
        mu.Wobs[(t-1)*R+i,j] <- 1/5*(W[(t-1)*R+i,j]+W[(t-0)*R+i,j]+W[(t+1)*R+i,j]+W[(t+2)*R+i,j]+W[(t+3)*R+i,j])
        W5obs[(t-1)*R+i,j] ~  dnorm(mu.Wobs[(t-1)*R+i,j], sd=W5se[(t-1)*R+i,j])
        
      }
      
    }
    
  }  
  
  # model for 1-year county estimates
  for(j in 1:gp) {
    
    for(i in 1:n0) {
      
      W1obs[i, j] ~ dnorm(W[ony_ind[i], j], sd=W1se[i, j])
      
    }
    
  }
  
  # model for state level 1-year estimates, mu.W defined from 2012-2021
  for(t in 1:n1){
    
    for(j in 1:gp){
      
      W1.state[t,j] ~ dnorm(mu.W[ony_ind_state[t],j], sd=W1.state.se[t,j])
      
    }
    
  }
  
  # model for W (W is a T*R by gp matrix, in this case 10*100 by 3 matrix, t starts in 2012 and ends in 2021)
  for(j in 1:gp){
    
    for(t in 1:(T+4)){ # T+4==2021 for W
      
      for(i in 1:R){
        
        W[(t-1)*R+i,j] ~ dnorm(mu.W[t,j] + g[(t-1)*R+i,j], tau=tau.W[j])
        W_ilogit[(t-1)*R+i,j] <- ilogit(W[(t-1)*R+i,j])
        
      }
      
      mu.W[t, j] ~ dnorm(mu.W0[j], tau.W0[j])
      
      # make g ICAR for each t, j; spatial random effect term for statewide average of variable j in year t
      g[((t-1)*R+1):((t-1)*R+R), j] ~ dcar_normal(adj=adj[], 
                                                  num=num[],
                                                  tau=tau.g[j],
                                                  zero_mean=1)
      
      
    }
    
    tau.g[j] ~ dgamma(.5, .5)
    tau.W[j] ~ dgamma(.5, .5)
    tau.W0[j] ~ dgamma(.5, .5)
    mu.W0[j] ~ dnorm(0, sd = 10)
    
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
  
  
  # fixed effects for data level
  for(p in 1:dp){
    
    beta.D[p] ~ dflat()
    
  }
  
  for(p in 1:tp){
    
    beta.T[p] ~ dflat()
    
  }
  
  
  for(p in 1:bp){
    
    beta.B[p] ~ dflat()
    
  }
  
  for(j in 1:gp) {
    
    gamma[j] ~ dflat() # fixed effects for process level ACS variables
    
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
yfit <- yfit # outcome matrix
n <- length(num) # number of NC counties
n0 <- nrow(logit_W1.obs) # number of counties/years with data in W1.obs
n1 <- nrow(logit_W1.state) # number of counties/years with data in W1.state
T <- length(2016:2021) # number of years
K <-  dim(yfit[, 1:3])[2] # number of outcomes

# prepare constraint objects
cens_index <- which(is.na(yfit$`People served by treatment programs`))
censored <- rep(1, 600)
yinit <- yfit[,1:3]
yinit[cens_index, 2] <-  2
c <- matrix(0, nrow=600, ncol=2)
c[cens_index, 1] <- 0
c[cens_index, 2] <- 4
c[-cens_index, 1] <- -Inf
c[-cens_index, 2] <- Inf

mod_constants <- list(R = n,
                      T = T,
                      K = K,
                      n0 = n0,
                      n1 = n1,
                      dp = ncol(X.D), # illicit opioid overdose deaths
                      tp = ncol(X.T), # people served by treatment programs
                      bp = ncol(X.B), # people receiving buprenorphine
                      gp = ncol(logit_W1.obs), # process level fixed effects (# of fixed effects for log(relative risk))
                      L = length(ell.rate),
                      ell.rate = ell.rate,
                      P = yfit$pop,
                      num = num,
                      adj = adj, 
                      S.se = logit_S.se,
                      mean.mu = rep(0, 2),
                      cov.mu = 10^4*diag(2),
                      mean.eps = rep(0, 3),
                      cov.eps.R = diag(3),
                      W5se = logit_W5.se,
                      W1se = logit_W1.se,
                      W1.state.se = logit_W1.state.se,
                      ony_ind = ony_ind,
                      ony_ind_state = ony_ind_state,
                      c = c
)

mod_data <- list(y=as.matrix(yfit[,1:3]),
                 censored=censored,
                 S=logit_S,
                 X.D=X.D,
                 X.T=X.T,
                 X.B=X.B,
                 W5obs = logit_W5.obs,
                 W1obs = logit_W1.obs,
                 W1.state = logit_W1.state
)

# specify initial values
beta.mu.init <- lm(S~ell.rate)$coefficients
logit_beta.mu.init <- lm(logit_S~ell.rate)$coefficients # use for beta.mu init 
beta.init <- matrix(data = 0, nrow = T, ncol = K)
Ninit <- floor(yfit$pop*beta.mu.init[1])
finit <- matrix(data = 0, nrow = n*T, ncol = K)
uinit <- rep(0, n*T)
ginit <- matrix(0, nrow = n*(T+4), ncol = ncol(logit_W1.obs))

# constraints on domain of various outcomes. Make sure inits for N are less than P and greater than the outcome counts
II <- which(Ninit > yfit$pop)

if(length(II) > 0){
  
  Ninit[II] <- floor(.1*yfit[II, 7])
  
}

II <- which(Ninit < yfit[,1])

if(length(II) > 0){
  
  Ninit[II] <- yfit[II,1] + 100
  
}

II <- which(Ninit < yfit[,2])

if(length(II) > 0){
  
  Ninit[II] <- yfit[II,2] + 100
  
}

II <- which(Ninit < yfit[,3])

if(length(II) > 0){
  
  Ninit[II] <- yfit[II,3] + 100
  
}

# set initial values.
mod_inits <- list(y = as.matrix(yinit),
                  N = Ninit,
                  beta = beta.init,
                  beta.D = rep(0,ncol(X.D)),
                  beta.T = rep(0,ncol(X.T)),
                  beta.B = rep(0,ncol(X.B)),
                  cov.eps = diag(3),
                  gamma = rep(0, ncol(logit_W1.obs)),
                  tau.v = .1, 
                  u = uinit,
                  f = finit,
                  g = ginit,
                  v = rep(0, n*T),
                  tau = rep(.1, K), 
                  phi.u = .5, 
                  phi.f = rep(.5, K), 
                  beta.mu = logit_beta.mu.init,
                  W = as.matrix(logit_W_init)
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
                                      'tau.g',
                                      'tau.W',
                                      'tau.W0',
                                      'eps',
                                      'cov.eps',
                                      'v',
                                      'beta',
                                      'beta.D',
                                      'beta.T',
                                      'beta.B',
                                      'gamma',
                                      'beta.mu',
                                      'N',
                                      'u',
                                      'f',
                                      'g',
                                      'pi',
                                      'lambda',
                                      'W',
                                      'W_ilogit',
                                      'mu',
                                      'mu.W',
                                      'mu.W0'
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

save(samples, file = "../output/full/NC_abundance_output_full_revision.Rda")

q()