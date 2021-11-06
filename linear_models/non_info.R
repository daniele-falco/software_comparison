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

#####################MODEL: linear model, non informative prior



####################SIMULATION OF DATA
set.seed(431234)
N=100 #number of observations
K<-4  #length of our parameters



# Matrix X of covariates with continuous covariates
X <- matrix(nrow=N, ncol=K)
X[,1] <- rep(1,N) #intercept
for(i in 2:K){
  X[,i] <- rnorm(n = N, mean = 0 ,sd = 1) #covariates 
}

# TRUE values
beta <- runif(K, min = -7, max = 7)
sigmasquare <- runif(1, min = 2, max = 10)
tau <- sigmasquare^(-1) 

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
#####################JAGS
model_NC<-
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
  for(j in 1:K){
    beta[j] ~ dnorm(0,0.0001) # stochastic node
    # remember that JAGS wants the precision
  }  
  sigma ~ dunif(0,1000) # stochastic node
  tau <- pow(sigma,-2) # deterministic node   
}"
model2<-textConnection(model_NC)

## Compilation phase
t<- system.time(
  {
    JAGS2<-jags.model(file=model2,
                      data=list('N'=N,'K'= p, 'X'=X,'Y'=Y),
                      inits=list('beta'=b0,'sigma'= sigma2),
                      n.chains=1,
                      n.adapt=5000)
    
    ## Sampling phase
    ## Run with CODA
    CODA2<-coda.samples(model=JAGS2,
                        variable.names=c('beta','tau', 'sigma'), 
                        n.iter=10000,
                        thin=2)
  })



#####################STAN

data_lm_nc <-list(N = N, 
                  y = Y,
                  p = p, 
                  X = as.matrix(X)
) 

# Initial values
inits_nc <- function() 
{
  list(
    beta = b0, 
    sigma=sigma2)
}

t<- system.time(
  {
    LM_NC <- stan(file = "LinearModel_non_informative.stan", 
                  data = data_lm_nc,
                  chains = 1, 
                  iter = 15000, 
                  warmup = 5000, 
                  thin= 2, 
                  seed = 42, 
                  init = inits_nc,
                  algorithm = 'NUTS')
  })


#####################NIMBLE

linearCode <- nimbleCode({
  for(i in 1:n){
    for(j in 1:k){
      a[i,j] <- X[i,j]*beta[j]
    }
    mu[i] <- sum(a[i,])
    # likelihood
    Y[i] ~ dnorm(mu[i],sd=sigma)
  }
  ## priors
  for(j in 1:k){
    beta[j] ~ dnorm(0,sd=100)
    
  }
  
  sigma ~ dunif(0.01,1000)
  
})

linearConsts <- list(n=N,k=p)
linearInits <-  list(beta=b0,sigma=sigma2)
linearData <- list(X=X,Y=Y)

t1<- system.time(
  {
    linear <- nimbleModel(code = linearCode,dimensions = list(a=c(N,p)),  
                          name = "linear",constants = linearConsts, data = linearData, inits = linearInits)
    Clinear <- compileNimble(linear)
    linearConf <- configureMCMC(linear, print = TRUE)
    linearConf$addMonitors(c("sigma", "beta"))
    
    linearMCMC <- buildMCMC(linearConf)
    ClinearMCMC <- compileNimble(linearMCMC, project = linear)
    t2<- system.time(
      {
        Psamples <- runMCMC(ClinearMCMC, niter=15000, nburnin=5000, thin=2, nchains=1,samplesAsCodaMCMC = TRUE,summary=TRUE)
      })
    
  })
