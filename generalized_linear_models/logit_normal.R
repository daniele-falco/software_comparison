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
#MODEL: Logistic regression, normal prior



#SIMULATION OF DATA
set.seed(431234)
N=1000 #number of observations
K<-16  #length of our parameters 16 
tm=NA
# Matrix X of covariates
X <- matrix(nrow=N, ncol=K)
X[,1] <- rep(1,N) #intercept
for(i in 2:K){
  X[,i] <- rnorm(n = N, mean = 0 ,sd = 1) #standardized covariates
}


beta <- runif(K, min = -7, max = 7)


# Y: response
Y <- rep(0,N)
for (i in 1:N){
  Y[i] <- rbern(1,1/(1+exp(-(X[i,]%*%beta))))
}

# Data
N <- length(Y)
p <- dim(X)[2]

# Initial values
b0 <- rep(0,K)
sigma0 <- 1 
#JAGS

# Compilation Model phase
model_C<-
  "model{
  for(i in 1:N){
    for (j in 1:K){
    a[i,j] <- X[i,j]*beta[j] 
    }
    mu[i]<-sum(a[i,]) # deterministic node
    ## likelihood
    Y[i] ~ dbern(1/(1+exp(-(mu[i]))))  # stochastic node
  }

  ## priors
  #REMEMBER Jags uses precision matrix
  
   for(j in 1:K){
    beta[j] ~ dnorm(0,0.1) # stochastic node
    
  } 
 
}"
model0<-textConnection(model_C)

## Compilation phase

t<- system.time(
  {
    JAGS0<-jags.model(file=model0,
                      data=list('N'=N,'K'= p, 'X'=X,'Y'=Y),
                      inits=list('beta'= beta), 
                      n.chains=1,
                      n.adapt=10000) 
    ## Sampling phase
    ## Run with CODA
    CODA0<-coda.samples(model=JAGS0,
                        variable.names=c('beta'), 
                        n.iter=10000,
                        thin=2)
  })


#STAN


inits <- function() 
{
  list(
    beta = b0)
} 

data_probit <-list(N = N, 
                   Y=Y,
                   p = p, 
                   X=as.matrix(X)
                   
                   
) 

t<- system.time(
  {
    PM <- stan(file = "logit_normal.stan", 
               data = data_probit,
               chains = 1, 
               iter = 20000, 
               warmup = 10000, 
               thin= 2, 
               seed = 42, 
               init = inits,
               algorithm = 'NUTS')
  })



#NIMBLE


probitCode <- nimbleCode({
  for(i in 1:N){
    for(j in 1:k){
      a[i,j] <- X[i,j]*beta[j]
    }
    mu[i] <- sum(a[i,])
    # likelihood
    Y[i] ~ dbern(1/(1+exp(-(mu[i]))))
  }
  ## priors
  for(i in 1:k){    
    b0[i] <- 0
  }
  
  for(i in 1:k){
    beta[i] ~ dnorm(0,sd=sqrt(10))
    
  }   
  
})

linearConsts <- list(N =N,k=p)
linearInits <-  list(beta=b0)
linearData <- list(X=X,Y=Y)


t1<- system.time(
  {
    linear <- nimbleModel(code = probitCode,dimensions = list(a=c(N,p)),
                          name = "linear",constants = linearConsts, data = linearData, inits = linearInits)
    Clinear <- compileNimble(linear)
    linearConf <- configureMCMC(linear, print = TRUE)
    linearConf$addMonitors(c("beta"))
    linearMCMC <- buildMCMC(linearConf)
    ClinearMCMC <- compileNimble(linearMCMC, project = linear)
    
    t2<- system.time(
      {
        Psamples <- runMCMC(ClinearMCMC, niter=20000, nburnin=10000, thin=2,
                            nchains=1,samplesAsCodaMCMC = TRUE,summary=TRUE)
      })
  })


