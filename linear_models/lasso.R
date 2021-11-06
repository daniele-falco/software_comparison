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

#MODEL: linear model, lasso prior



#SIMULATION OF DATA
set.seed(431234)
N=100 #number of observations
K<-16  #length of our parameters 16 or also 120 

# Matrix X of covariates
X <- matrix(nrow=N, ncol=K)
X[,1] <- rep(1,N) #intercept
for(i in 2:K){
  X[,i] <- rnorm(n = N, mean = 0 ,sd = 1) #covariates 
}

#when p=16
beta=c(1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5,0.01)

#when p=30
#beta=c(1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5,
#       0.01,1.2,-1.5,-2.9,3.2,1.4,0.05,0.5,-0.009,0.001,2.5,4.5,-0.09,0,0)
#beta=c(1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#beta=c( 2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#when p=100
#beta=c(1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5,
#       0.01,1.2,-1.5,-2.9,3.2,1.4,0.05,0.5,-0.009,0.001,2.5,4.5,-0.09,1,2,3,-1,-2,-3,0.5,1.5,2.5,-1,-1.5,-2.5,-3,0.2,0.6,0.7,-0.5,-0.9,1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5,
#       0.01,1.2,-1.5,-2.9,3.2,1.4,0.05,0.5,-0.009,0.001,2.5,4.5,-0.09,1,2,3,-1,-2,-3,0.5,1.5,2.5,-1,-1.5,-2.5,-3,0.2,0.6,0.7,-0.5,-0.9,-0.2,0.3,0.9,1.5,1.8,2.4,0,0)
#beta=c(1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5,
#       0.01,1.2,-1.5,-2.9,3.2,1.4,0.05,0.5,-0.009,0.001,2.5,4.5,-0.09,1,2,3,-1,-2,-3,0.5,1.5,2.5,-1,-1.5,-2.5,-3,0.2,0.6,0.7,-0.5,-0.9,1,-1,-2,3,
#       0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0)
#beta=c( 2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#for repeated simulations with p=30
#beta <- runif(15, min = -7, max = 7)
#beta=c(beta,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)


sigmasquare <- 5
# for repeated simulations
#sigmasquare <- runif(1, min = 2, max = 10)
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
sigma2 <- 1
lambda2=1
sigma0 <- 1 
nu0 <- 0.0001 

#JAGS
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
  for (j in 1:K){
   beta[j] ~ ddexp(0.0,  (pow(lambda2,0.5))) #attenzione, stan e jags usano parametrizzazioni diverse
  }
  lambda2 ~ dexp(0.1)
  tau ~ dgamma(nu0/2,nu0*sigma0*sigma0/2) # stochastic node
  
}"
model0<-textConnection(model_C)

## Compilation phase
t<- system.time(
  {
    JAGS0<-jags.model(file=model0,
                      data=list('N'=N,'K'= p, 'X'=X,'Y'=Y,'nu0'= nu0, 'sigma0' = sigma0, 'lambda2'=lambda2),
                      inits=list('beta'= beta,'tau'= 1/sigma2), 
                      n.chains=1,
                      n.adapt=10000) 
    ## Sampling phase
    ## Run with CODA
    CODA0<-coda.samples(model=JAGS0,
                        variable.names=c('beta','tau', 'sigmasq'), 
                        n.iter=10000,
                        thin=2)
  })


#STAN

inits <- function() 
{
  list(
    beta = b0, 
    sigma2=sigma2,
    lambda2=lambda2)
} 

data_lm_zellner <-list(N = N, 
                       y=Y,
                       p = p, 
                       X=as.matrix(X),
                       
                       b0=rep(0,p),
                       sigma0=sigma0,
                       nu0=nu0
) 

t<- system.time(
  {
    LM_Z1 <- stan(file = "lasso.stan", 
                  data = data_lm_zellner,
                  chains = 1, 
                  iter = 20000, 
                  warmup = 10000, 
                  thin= 2, 
                  seed = 42, 
                  init = inits,
                  algorithm = 'NUTS')
  })


#NIMBLE


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
  for (j in 1:k){
    beta[j] ~ ddexp(location=0, scale=1.0 / (pow(lambda2,0.5)))
  }
  lambda2 ~ dexp(0.1)
})

linearConsts <- list(n =N,k=p)
linearInits <-  list(beta=b0,sigmasq=sigma2,tau = 1/sigma2, lambda2=lambda2)
linearData <- list(X=X,Y=Y)


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
        Psamples <- runMCMC(ClinearMCMC, niter=20000, nburnin=10000, thin=2,
                            nchains=1,samplesAsCodaMCMC = TRUE,summary=TRUE)
      })
  })