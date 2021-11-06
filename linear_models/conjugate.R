library(coda)
library(tictoc)
library(mvtnorm)

# Jags
library(rjags)

# Stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE))

# Nimble 
library(nimble) 
library(ggplot2)
library(extraDistr)
#####################MODEL: linear model, conjugate prior



####################SIMULATION OF DATA
set.seed(1)
N=100 #number of observations
K<-4  #length of our parameters

# Matrix X of covariates with binary covariates
X <- matrix(nrow=N, ncol=K)
X[,1] <- rep(1,N) #intercept

if(K==4){
  p_b=c(0.1,0.5,0.8)}
if(K==16){
  p_b=c(0.1,0.2,0.3,0.4,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.7,0.8,0.9)}
if(K==50){
  p_b=runif((K-1),min=0.05,max=0.95)}

for(i in 2:K){
  X[,i] <- rbern(N,p_b[i-1]) #covariates 
}

# Matrix X of covariates with continuous covariates
X <- matrix(nrow=N, ncol=K)
X[,1] <- rep(1,N) #intercept
for(i in 2:K){
  X[,i] <- rnorm(n = N, mean = 0 ,sd = 1) #covariates 
}

# TRUE values
beta <- runif(K, min = -7, max = 7)
sigmasquare <- runif(1, min = 2, max = 10)
tau <- sigmasquare^(-1) #Since Jags uses precision rather than variance

# Y: response
Y <- rep(0,N)
for (i in 1:N){
  Y[i] <- rnorm(n=1, mean=X[i,]%*%beta,sd=sigmasquare^0.5)
}


# Data
N <- length(Y)
p <- dim(X)[2]

# Initial values
b0 <- rep(0,K)
sigma0 <- 1 
nu0 <- 0.0001 
B0 <- diag(p)

#####################JAGS
model_C<-
  "model{
  for(i in 1:N){
    for (j in 1:K){
    a[i,j] <- X[i,j]*beta[j] 
    }
    mu[i]<-sum(a[i,]) # deterministic node
    ## likelihood
    Y[i] ~ dnorm(mu[i],tau)  # stochastic node
  }

  ## priors
  #REMEMBER Jags uses precision matrix
  sigmasq <- pow(tau,-1) # deterministic node
  
  beta[1:K] ~ dmnorm(rep(0,K),inverse(sigmasq*B0[1:K,1:K])) # stochastic node
  tau ~ dgamma(nu0/2,nu0*(sigma0*sigma0)/2) # stochastic node
  
}"
model0<-textConnection(model_C)

## Compilation phase
t<- system.time(
  {
    JAGS0<-jags.model(file=model0,
                      data=list('N'=N,'K'= p, 'X'=X,'Y'=Y,'nu0'= nu0, 'sigma0' = sigma0, 'B0' = B0),
                      inits=list('beta'= b0,'tau'= 1/sigma0), 
                      n.chains=1,
                      n.adapt=1000) 
    ## Sampling phase
    ## Run with CODA
    CODA0<-coda.samples(model=JAGS0,
                        variable.names=c('beta','tau', 'sigmasq'), 
                        n.iter=10000,
                        thin=2)
  })


#####################STAN
inits <- function() 
{
  list(
    beta = b0, 
    sigma2=sigma0)
} 

data_lm_zellner <-list(N = N, 
                       y=Y,
                       p = p, 
                       X=as.matrix(X),
                       B0=B0,
                       b0=rep(0,p),
                       sigma0=sigma0,
                       nu0=nu0
) 

t<- system.time(
  {
    LM_Z1 <- stan(file = "LinearModel_conjugate.stan", 
                  data = data_lm_zellner,
                  chains = 1, 
                  iter = 11000, 
                  warmup = 1000, 
                  thin= 2, 
                  seed = 42, 
                  init = inits,
                  algorithm = 'NUTS')
  })


#####################NIMBLE

linearCodeZ <- nimbleCode({
  for(i in 1:n){
    for(j in 1:k){
      a[i,j] <- X[i,j]*beta[j]
    }
    mu[i] <- sum(a[i,])
    # likelihood
    Y[i] ~ dnorm(mu[i],sd=sigma)
  }
  ## priors
  for(i in 1:k){    
    b0[i] <- 0
  }
  
  nu0 <- 0.001
  sigma0 <- 1
  
  sigmasq ~ dinvgamma(nu0*0.5,nu0*sigma0*sigma0*0.5)
  tau <- pow(sigmasq,-1)  ## residual std dev
  sigma <- pow(sigmasq,0.5)
  cova[1:k,1:k] <- sigmasq*B0[1:k,1:k]
  beta[1:k] ~ dmnorm(b0[1:k],cov=cova[1:k,1:k])   
  
})

linearConsts <- list(n =N,k=p)
linearInits <-  list(beta=b0,sigmasq=sigma0,tau = 1/sigma0)
linearData <- list(X=X,Y=Y,B0=B0)


t1<- system.time(
  {
    linear <- nimbleModel(code = linearCodeZ,dimensions = list(a=c(N,p)),
                          name = "linear",constants = linearConsts, data = linearData, inits = linearInits)
    Clinear <- compileNimble(linear)
    linearConf <- configureMCMC(linear, print = TRUE)
    linearConf$addMonitors(c("tau", "sigmasq", "beta"))
    linearMCMC <- buildMCMC(linearConf)
    ClinearMCMC <- compileNimble(linearMCMC, project = linear)
    t2<- system.time(
      {
        Psamples <- runMCMC(ClinearMCMC, niter=11000, nburnin=1000, thin=2,
                            nchains=1,samplesAsCodaMCMC = TRUE,summary=TRUE)
      })
  })
