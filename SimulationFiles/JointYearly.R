
###################################################################################################################################
##################### Implement MCMC algorithm for the proposed joint model with yearly survey data cases #########################
####################### Code runs in a 4-iteration loop with 25 data sets in parallel in each iteration ###########################
###################################################################################################################################

library(iterators)
library(doParallel)
library(foreach)
library(doRNG)

for(MM in 1:4){
  filename = paste("JointYearly.", toString(MM), ".Rda", sep="")
  outfile = paste("JointYearlyLog.",toString(MM),".txt",sep="")
  Mrep = ((MM-1)*25+1):(MM*25) ### identify the set of 25 simulated data sets that are run in this code
  
  cl <- makePSOCKcluster(25) ## this code is designed to use 25 cores of a cluster
  registerDoParallel(cl)
  
  # Setup -------------------------------------------------------------------
  
  writeLines(c(""), outfile)
  set.seed(1025)
  
  SimMCMC<-foreach(kk=1:25) %dorng% {
    
    k=Mrep[kk]
    library(nimble) 
    library(coda) 
    
    sink(outfile, append=TRUE)
    cat(paste("Starting iteration",k,"\n"))
    sink()
    
    ############################################################ 
    ###########             Load data               ############
    ############################################################ 
    
    ### time intervals 
    TT=10
    n=100
    
    ### load in data
    load("SimulatedData.Rda")
    Ntrue = N[1:(n*TT),k]
    ftrue = f[1:(n*TT),k]
    gtrue = g[1:(n*TT),k]
    utrue = u[1:(n*TT),k]
    
    Y=cbind(D[1:(TT*n),k],T[1:(TT*n),k]) ### observed outcome counts
    S0 = S[1:(TT),k] ### survey data
    P=P[1:(TT*n)] ### population
    T=TT ## number of years
    
    ############################################################
    ###########           Set up            ############
    ############################################################
    ### information on neighborhood structure needed for CAR model
    num<-colSums(A)
    adj<-NULL
    for(j in 1:n){
      adj<-c(adj,which(A[j,]==1))
    }
    adj<-as.vector(adj)
    num<-as.vector(num)
    weights<-1+0*num
    n_regions=length(num) 
    
    ###############################################################################
    ######### Write the model #################
    ########################################################
    model_code=nimbleCode({
      for(t in 1:1){
        for(i in 1:R){
          pi.D[(t-1)*R+i] <- ilogit(inprod(X[(t-1)*R+i,],beta.D[])+(f[(t-1)*R+i]+mu.f[(t-1)*R+i])+eps.D[(t-1)*R+i])
          pi.T[(t-1)*R+i] <- ilogit(inprod(X[(t-1)*R+i,],beta.T[])+(g[(t-1)*R+i]+mu.g[(t-1)*R+i])+eps.T[(t-1)*R+i])
          mu.u[(t-1)*R+i] <- inprod(W[(t-1)*R+i,],gamma[])
          lambda[(t-1)*R+i] <- exp((u[(t-1)*R+i]+mu.u[(t-1)*R+i])+v[(t-1)*R+i])
          y1[(t-1)*R+i] ~ dbinom(pi.D[(t-1)*R+i],N[(t-1)*R+i])
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
        S[t] ~ T(dnorm(mu[t],sd=S.se[t]),0,1)
      }
      for(t in 2:T){
        for(i in 1:R){
          pi.D[(t-1)*R+i] <- ilogit(inprod(X[(t-1)*R+i,],beta.D[])+(f[(t-1)*R+i]+mu.f[(t-1)*R+i])+eps.D[(t-1)*R+i])
          pi.T[(t-1)*R+i] <- ilogit(inprod(X[(t-1)*R+i,],beta.T[])+(g[(t-1)*R+i]+mu.g[(t-1)*R+i])+eps.T[(t-1)*R+i])
          lambda[(t-1)*R+i] <- exp((u[(t-1)*R+i]+mu.u[(t-1)*R+i])+v[(t-1)*R+i])
          y1[(t-1)*R+i] ~ dbinom(pi.D[(t-1)*R+i],N[(t-1)*R+i])
          y2[(t-1)*R+i] ~ dbinom(pi.T[(t-1)*R+i],N[(t-1)*R+i])
          N[(t-1)*R+i] ~ dbinom(mu[t]*lambda[(t-1)*R+i],P[(t-1)*R+i])
          v[(t-1)*R+i] ~ dnorm(0,tau=tau.v)
          eps.D[(t-1)*R+i] ~ dnorm(0,tau=tau.D)
          eps.T[(t-1)*R+i] ~ dnorm(0,tau=tau.T)
          mu.f[(t-1)*R+i] <- phi.f*f[(t-2)*R+i]
          mu.g[(t-1)*R+i] <- phi.g*g[(t-2)*R+i]
          mu.u[(t-1)*R+i] <- inprod(W[(t-1)*R+i,],gamma[])+phi.u*(u[(t-2)*R+i]-inprod(W[(t-2)*R+i,],gamma[]))
        }
        mu[t] <- beta[1]+beta[2]*t
        f[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.f,zero_mean=1)
        u[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.u,zero_mean=1)
        g[((t-1)*R+1):((t-1)*R+R)] ~ dcar_normal(adj=adj[], num=num[], tau=tau.g,zero_mean=1)
        S[t] ~ T(dnorm(mu[t],sd=S.se[t]),0,1)
      }
      for(i in 1:bp){
        beta.T[i] ~ dflat()
        beta.D[i] ~ dflat()
      }
      for(i in 1:2){
        beta[i] ~ dflat()
        gamma[i] ~ dflat()
      }
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
    
    
    mod_constants=list(bp = length(X[1,]),R=n,T=TT,P=P,adj=adj,X=X,num=num,W=W,S.se=rep(.005,TT))
    mod_data=list(y1=Y[,1],y2=Y[,2],S=S0)
    
    ##### Set initial values for Markov chains
    
    Ninit = floor(1/2*(5*Y[,2]+10*Y[,1])+100) ### check that N is not greater than P and also not less than T or D
    II=which(Ninit>P)
    if(length(II)>0){
      Ninit[II]=floor(.1*P[II])
    }
    II = which(Ninit<Y[,1])
    if(length(II)>0){
      Ninit[II] = Y[II,1]+100
    }
    II = which(Ninit<Y[,2])
    if(length(II)>0){
      Ninit[II] = Y[II,2]+100
    }
    
    mod_inits=list(N=Ninit,beta.T=rep(-3,T),beta.D=rep(-6,T),beta=c(.05,-.001),tau.v=.1,tau.f=.5,tau.g=.5,tau.u=.5,
                   u=rep(0,n*T),v=rep(0,n*T),f=rep(0,n*T),g=rep(0,n*T),tau.D=.1,tau.T=.1,eps.T=rep(0,n*T),eps.D=rep(0,n*T),
                   gamma=rep(0,length(W[1,])),phi.u=.5,phi.f=.5,phi.g=.5)
    
    # Build the model.
    nimble_model <- nimbleModel(model_code, mod_constants,mod_data,mod_inits)
    compiled_model <- compileNimble(nimble_model,resetFunctions = TRUE)
    
    # Set up samplers.
    mcmc_conf <- configureMCMC(nimble_model,monitors=c('beta.T','beta.D','beta','N','u','f','g','gamma','lambda'),control=list(adaptive=TRUE,adaptInterval=100,sliceWidths=rep(100,T),sliceMaxSteps=10000,maxContractions=10000000,sliceAdaptWidthMaxIter=0,sliceAdaptFactorMaxIter=0),useConjugacy = TRUE)
    
    #### use AF_Slice update for N_{i,1},...,N_{i,T}
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
    samples=runMCMC(compiled_mcmc,inits=mod_inits,
                    nchains = 1, nburnin=MCS/2,niter = MCS,samplesAsCodaMCMC = TRUE,thin=50,
                    summary = FALSE, WAIC = FALSE,progressBar=TRUE) 
    
    ### organize output
    list(colMeans(samples),apply(samples,2,quantile,probs=c(.025,.975)),apply(samples,2,sd))
    
  }
  
  save(SimMCMC,file=filename)
}
stopCluster(cl)
