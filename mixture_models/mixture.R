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

#MODEL: MIXTURE MODEL



#SIMULATION OF DATA
set.seed(431234)
N=1000 #number of observations

H=4
#H=2



if(H==2){
  mu=mu=c(runif(1, min = -2, max = 0),runif(1, min = 1, max = 3))
  Y=c(rnorm(N/2,mean=mu[1],sd=1),rnorm(N/2,mean=mu[2],sd=1))
}


if(H==4){
  mu=c(-4,0,2,6)
  Y=c(rnorm(N/4,mean=mu[1],sd=1),rnorm(N/4,mean=mu[2],sd=1),rnorm(N/4,mean=mu[3],sd=1),rnorm(N/4,mean=mu[4],sd=1))
  
}


sigmasquare <- 1
tau <- sigmasquare^(-1) 


N <- length(Y)
# Initial values
sigma2 <- 1



#JAGS


model_C = "
model {
for (i in 1:N) {
Y[i] ~ dnorm(mu[zeta[i]], tau[zeta[i]])
zeta[i] ~ dcat(pi[])
}
v2~dgamma(1,1)
for (i in 1:H) {
mu0[i] ~ dnorm(0,v2)#precision di v2
tau[i] ~ dgamma(1,1)

sigma[i] <- 1/sqrt(tau[i])
}
mu[1:H] <- sort(mu0)
pi ~ ddirich(a)
}"
model0<-textConnection(model_C)

## Compilation phase
t<- system.time(
  {
    JAGS0<-jags.model(file=model0,
                      data=list('N'=N,'Y'=Y,'H'=H,'a'=rep(1,H)),
                      inits=list('tau'= rep(1/sigma2,H),'v2'=1),
                      n.chains=1,
                      n.adapt=10000) 
    ## Sampling phase
    ## Run with CODA
    CODA0<-coda.samples(model=JAGS0,
                        variable.names=c('mu','tau','v2'), 
                        n.iter=10000,
                        thin=2)
  })

#STAN



inits <- function() 
{
  list(
    theta = rep(1/H,H), 
    sigma2=rep(1,H),
    mu=seq(1,H,1),v=1)
} 

data_lm_zellner <-list(N = N, 
                       Y=Y,
                       H = H
) 

t<- system.time(
  {
    LM_Z1 <- stan(file = "mixture_model.stan", 
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

nimbleOptions(allowDynamicIndexing = TRUE)
linearCodeZ <- nimbleCode({
  
  p[1:H]~ddirich(a[1:H])
  for(i in 1:n) {
    z[i]~ dcat(p[1:H])
    y[i] ~ dnorm(mu[z[i]], tau=tau[z[i]])
    
  }
  for (h in 1:H) {
    mu[h]~ dnorm(0, var=v)
    tau[h]~dgamma(shape=1,rate=1)
    #sigma[i] <- 1/sqrt(tau[i])
  }
  
  
  v~dinvgamma(shape=1,rate=1)
  
  
})

linearConsts <- list(n =N, H=H
                     , a=rep(1,H)
)
linearInits <- list(mu=seq(1,H,1),
                    tau=rep(1,H),
                    z = rep(1,N)
                    
)
linearData <- list(y=Y)


t1<- system.time(
  {
    rModel <- nimbleModel(linearCodeZ,dimensions = list(p = H, a = H, mu=H, tau=H), data = linearData, inits = linearInits, constants = linearConsts)
    
    cModel <- compileNimble(rModel)
    conf <- configureMCMC(rModel, monitors = c("mu", "tau","v"))
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project = rModel)
    t2<- system.time(
      {
        Psamples <- runMCMC(cmcmc, niter=20000, nburnin=10000, thin=2,
                            nchains=1,samplesAsCodaMCMC = TRUE,summary=TRUE)
      })
  })
