############## Simulate Data for Simulation Study of R21
############## Assume we have survey data at each year
library(ngspatial)
library(MASS)



# Run Code ----------------------------------------------------------------


M=100 ### Number of simulation data sets
TT=10 ### Number of years
n=100 ### number of counties (10x10 grid)
L = TT ### amount of survey information
ell = seq(1,TT,by=1) ### the years with survey information

A = adjacency.matrix(sqrt(n),sqrt(n)) ### neighborhood matrix
H = diag(rowSums(A))

### simulate the population size for each year
### create nT-dim vector with (year1, year2, ...,  yearT)
set.seed(113)
P = rep(0,n*TT)
t = 1
Indt = ((t-1)*n+1):(t*n)
mm0 = 120000
ss0 = 250000
lPop = rnorm(n,log(mm0)-1/2*log(ss0^2/(mm0^2)+1),sqrt(log(ss0^2/(mm0^2+1))))
P[Indt] = ceiling(exp(lPop))
for(t in 2:TT){
  Indt = ((t-1)*n+1):(t*n)
  lPop = rnorm(n,lPop,.1)
  P[Indt] = ceiling(exp(lPop))
}

### create design matrix X that has different intercepts for each year
X = matrix(0,n*TT,TT)
for(i in 1:TT){
  X[((i-1)*n+1):((i-1)*n+n),i]=rep(1,n)
}


### create matrix W with two spatially varying covariates
### must be centered at zero and no intercept
W = cbind(rep(0,n*TT),rep(0,n*TT))
t=1
Indt = ((t-1)*n+1):(t*n)
temp = mvrnorm(1,rep(0,n),.5*solve(H-.8*A)) ### covariate that slightly changes over time
temp2 = mvrnorm(1,rep(0,n),2*solve(H-.5*A)) ### covariate that slightly changes over time
W[Indt,1]=temp-mean(temp)
W[Indt,2]=temp2-mean(temp2)
for(t in 2:TT){
  Indt = ((t-1)*n+1):(t*n)
  Indt0 = ((t-2)*n+1):((t-1)*n)
  W[Indt,1] = W[Indt0,1]+rnorm(n,0,.25) ### covariate that slightly changes over time
  W[Indt,1] = W[Indt,1]-mean(W[Indt,1])
  W[Indt,2] = W[Indt0,2] + rnorm(n,0,.25)
  W[Indt,2] = W[Indt,2]-mean(W[Indt,2])
}


### set values of hyperparameters

sig.v = .005
sig.T = .001
sig.D = .005
phi.u = .7
phi.f = .5
phi.g = .7
phi.a = .4
tau.u = .2
tau.f = .15
tau.g = .2
beta = c(.05,-.0001)
beta.T = c(rep(-3,4),rep(-2.5,3),rep(-2,3))
beta.D = c(rep(-5,5),rep(-4.5,2),rep(-4,3))
gamma = c(.35,-.25)


### empty matrices to store simulated values
N = matrix(0,n*TT,M)
u = matrix(0,n*TT,M)
v = matrix(0,n*TT,M)
g = matrix(0,n*TT,M)
f = matrix(0,n*TT,M)
eps.T = matrix(0,n*TT,M)
eps.D = matrix(0,n*TT,M)
T = matrix(0,n*TT,M)
D = matrix(0,n*TT,M)
S = matrix(0,L,M)

set.seed(516219)
for(m in 1:M){
  
  ### simulate f 
  t=1
  Indt = ((t-1)*n+1):(t*n)
  f0 = mvrnorm(1,rep(0,n),tau.f*solve(H-.95*A)) ## simulate with proper CAR model
  f[Indt,m] = f0-mean(f0)
  for(t in 2:TT){
    Indt = ((t-1)*n+1):(t*n)
    Indt0 = ((t-2)*n+1):((t-1)*n)
    f0 = mvrnorm(1,phi.f*f[Indt0,m],tau.f*solve(H-.95*A)) ## simulate with proper CAR model
    f[Indt,m] = f0-mean(f0)
  }
  
  ### simulate g 
  t=1
  Indt = ((t-1)*n+1):(t*n)
  g0 = mvrnorm(1,rep(0,n),tau.g*solve(H-.95*A)) ## simulate with proper CAR model
  g[Indt,m] = g0-mean(g0)
  for(t in 2:TT){
    Indt = ((t-1)*n+1):(t*n)
    Indt0 = ((t-2)*n+1):((t-1)*n)
    g0 = mvrnorm(1,phi.g*g[Indt0,m],tau.g*solve(H-.95*A)) ## simulate with proper CAR model
    g[Indt,m] = g0-mean(g0)
  }
  
 ### simulate u
  t=1
  Indt = ((t-1)*n+1):(t*n)
  u0 = mvrnorm(1,W[Indt,]%*%gamma,tau.u*solve(H-.8*A)) ## simulate with proper CAR model
  u[Indt,m] = u0-mean(u0)
  for(t in 2:TT){
    Indt = ((t-1)*n+1):(t*n)
    Indt0 = ((t-2)*n+1):((t-1)*n)
    u0 = mvrnorm(1,W[Indt,]%*%gamma+phi.u*(u[Indt0,m]-W[Indt0,]%*%gamma),tau.u*solve(H-.8*A)) ## simulate with proper CAR model
    u[Indt,m] = u0-mean(u0)
  }
  
  for(t in 1:TT){
    Indt = ((t-1)*n+1):(t*n)
    v0 = rnorm(n,0,sqrt(sig.v))
    v[Indt,m] = v0
  }
  
  Z = u[,m]+v[,m]
  lambda = exp(Z)
  
  
  ### simulate N
  for(t in 1:TT){
    Indt = ((t-1)*n+1):(t*n)
    N[Indt,m] = rbinom(n,P[Indt],(beta[1]+beta[2]*t)*lambda[Indt])
  }
 
  
  ### simulate epsT and epsD
  for(t in 1:TT){
    Indt = ((t-1)*n+1):(t*n)
    eT0 = rnorm(n,0,sqrt(sig.T))
    eD0 = rnorm(n,0,sqrt(sig.D))
    eps.T[Indt,m] = eT0 
    eps.D[Indt,m] = eD0 
  }
  
  
  ### simulate survey data 
  s.se = .005
  for(ll in 1:L){
    l = ell[ll]
    S[l,m] = rnorm(1,(beta[1]+beta[2]*l),s.se)
  }

  
  ### simulate treatment and death data
  logitp = X%*%beta.D+f[,m]+eps.D[,m]
  logitq = X%*%beta.T+g[,m]+eps.T[,m]
  T[,m] = rbinom(n*TT,N[,m],exp(logitq)/(1+exp(logitq)))
  D[,m] = rbinom(n*TT,N[,m],exp(logitp)/(1+exp(logitp)))
  
  
}

save(X,W,A,P,T,D,S,u,v,N,f,g,eps.T,eps.D,file = "SimulatedData.Rda")



