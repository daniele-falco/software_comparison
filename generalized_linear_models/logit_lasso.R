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

#MODEL: Logistic regression, lasso



#SIMULATION OF DATA
set.seed(431234)
N=1000 #number of observations
K<-16  #length of our parameters 16 or also 120 

# Matrix X of covariates
X <- matrix(nrow=N, ncol=K)
X[,1] <- rep(1,N) #intercept
for(i in 2:K){
  X[,i] <- rnorm(n = N, mean = 0 ,sd = 1) #covariates 
}



# TRUE values
#when p=16
beta=c(1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5,0.01)

#when p=100
#beta=c(1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5,
#       0.01,1.2,-1.5,-2.9,3.2,1.4,0.05,0.5,-0.009,0.001,2.5,4.5,-0.09,1,2,3,-1,-2,-3,0.5,1.5,2.5,-1,-1.5,-2.5,-3,0.2,0.6,0.7,-0.5,-0.9,1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5,
#       0.01,1.2,-1.5,-2.9,3.2,1.4,0.05,0.5,-0.009,0.001,2.5,4.5,-0.09,1,2,3,-1,-2,-3,0.5,1.5,2.5,-1,-1.5,-2.5,-3,0.2,0.6,0.7,-0.5,-0.9,-0.2,0.3,0.9,1.5,1.8,2.4,0,0)
#beta=c(1,-1,-2,3,2,0.5,0.05,-0.005,0.001,1,5,-0.02,2.5,1.5,-0.5,
#       0.01,1.2,-1.5,-2.9,3.2,1.4,0.05,0.5,-0.009,0.001,2.5,4.5,-0.09,1,2,3,-1,-2,-3,0.5,1.5,2.5,-1,-1.5,-2.5,-3,0.2,0.6,0.7,-0.5,-0.9,1,-1,-2,3,
#       0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0)
#beta=c( 2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)



sigmasquare <- 5
tau <- sigmasquare^(-1) #Since Jags uses precision rather than variance


# Y: response
Y <- rep(0,N)
for (i in 1:N){
  Y[i] <- rbern(1,1/(1+exp(-(X[i,]%*%beta))))
}


N <- length(Y)
p <- dim(X)[2]

# Initial values
b0 <- rep(0,K)

#JAGS
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

  for (j in 1:K){
   beta[j] ~ ddexp(0.0,  (pow(lambda2,0.5)))
  }
  lambda2 ~ dexp(0.1)
 
}"
model0<-textConnection(model_C)

## Compilation phase
t<- system.time(
  {
    JAGS0<-jags.model(file=model0,
                      data=list('N'=N,'K'= p, 'X'=X,'Y'=Y),
                      inits=list('beta'= b0,'lambda2'=1), 
                      n.chains=1,
                      n.adapt=10000) 
    ## Sampling phase
    ## Run with CODA
    CODA0<-coda.samples(model=JAGS0,
                        variable.names=c('beta'), 
                        n.iter=5000,
                        thin=2)
  })


#STAN
inits <- function() 
{
  list(
    beta = b0,
    lambda2=1)
} 

data_probit <-list(N = N, 
                   Y=Y,
                   p = p, 
                   X=as.matrix(X)
                   
) 

t<- system.time(
  {
    PM <- stan(file = "logit_model_lasso.stan", 
               data = data_probit,
               chains = 1, 
               iter = 15000, 
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
  for(i in 1:k){    # 
    b0[i] <- 0
  }
  
  for (j in 1:k){
    beta[j] ~ ddexp(location=0, scale=1.0 / (pow(lambda2,0.5)))
  }
  lambda2 ~ dexp(0.1)
  
  
})

linearConsts <- list(N =N,k=p)
linearInits <-  list(beta=b0,lambda2=1)
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
        Psamples <- runMCMC(ClinearMCMC, niter=15000, nburnin=10000, thin=2,
                            nchains=1,samplesAsCodaMCMC = TRUE,summary=TRUE)
      })
  })


