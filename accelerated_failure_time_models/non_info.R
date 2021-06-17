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

#MODEL: AFT, non-informative prior



#SIMULATION OF DATA

N <- 1000 # number of elements in the dataset
# simulate covariates 

p<-16
X<-matrix(nrow=N, ncol=p)
X[,1]<-rep(1,N) # intercept
for(i in 2:p){
  X[,i] <- rnorm(n = N, mean = 0 , sd =1 ) 
}

# true parameters

beta <- runif(p,min=-1,max=1) 
mu <- X %*% beta

sigma <- 1
#for repeated simulations
sigma <- runif(1, min = 1, max = 5)

alpha <- 1/sigma 
lambda <- exp(mu)/((log(2))^(1/alpha)) #attenzione alle diverse parametrizzazione della weibull
perc=0.2 #Percentuale dati censurati. Oppure 0.5 o 0.8
#perc=runif(1, min = 0.2, max = 0.8) #simulazioni ripetute
theta=(1/perc-1)^(1/alpha)*lambda


t = rweibull(N, shape=alpha, scale = lambda) #t_i
cent = rweibull(N, shape=alpha, scale = theta) #c_i

censured=t>cent
delta <- as.logical(1-censured)
t[censured==1]=NA
cent[censured==0]=Inf

beta0=rep(0.1, p)
alpha0=0.1


#JAGS

model_AFT1<-"model {
 for(i in 1:N){
   for(s in 1:p){
   a[i,s] <- X[i,s]*(beta[s])
   }
   mu[i] <- sum(a[i,])
   lambda[i] <- log(2)*exp(-mu[i]*alpha) 
   # likelihood
   censured[i] ~ dinterval(t[i], cent[i])
   t[i] ~ dweib(alpha, lambda[i])
   }
   #Priors
   for(i in 1:p){
   beta[i] ~ dnorm(0,1/1000)
   }
   alpha ~ dunif(0.01,100)
}"
model1<-textConnection(model_AFT1)

## Compilation phase
t<- system.time(
  {
    JAGS1<-jags.model(file=model1,
                      data=list('N'=N,'X'=X,'p'= p,'cent'=cent, 't'=t,'censured'=censured),
                      inits=list('beta'=beta0,'alpha'=alpha0),
                      n.chains=1,
                      n.adapt=5000) 
    ## Sampling phase
    ## Run with CODA
    CODA1<-coda.samples(model=JAGS1,
                        variable.names=c('beta','alpha'), 
                        n.iter=5000,
                        thin=2) 
  })




#NIMBLE


AFT_NH_Code <- nimbleCode({
  for(i in 1:N){
    for(s in 1:p){
      a[i,s] <- X[i,s]*(beta[s])
    }
    mu[i] <- sum(a[i,])
    lambda[i] <- log(2)*exp(-mu[i]*alpha) 
    # likelihood
    censured[i] ~ dinterval(t[i], cent[i])
    t[i] ~ dweib(alpha, lambda[i])
  }
  #Priors
  for(i in 1:p){
    beta[i] ~ dnorm(0,var=1000)
  }
  alpha ~ dunif(0.01,100)
})



aftConsts <- list(N=N,p=p)
aftInits <-  list(beta=beta0,alpha=alpha0)
aftData <- list(X=X,cent=cent,t=t,censured=censured)

t1<- system.time(
  {
    aftNH <- nimbleModel(code = AFT_NH_Code,dimensions = list(a = c(N,p)),
                         name = "aftNH",constants = aftConsts, data = aftData, inits = aftInits)
    Caft <- compileNimble(aftNH)
    aftConf <- configureMCMC(aftNH, print = TRUE)
    aftConf$addMonitors(c("beta","alpha"))
    aftMCMC <- buildMCMC(aftConf)
    CaftMCMC <- compileNimble(aftMCMC, project = aftNH)
    t2<- system.time(
      {
        Psamples <- runMCMC(CaftMCMC, niter=(10000), nburnin=5000, thin=2,
                            nchains=1,samplesAsCodaMCMC = TRUE,summary=TRUE)
      })
  })








####SIMULATION OF DATA

N <- 1000 # number of elements in the dataset
# simulate covariates 

p<-16
X<-matrix(nrow=N, ncol=p)
X[,1]<-rep(1,N) # intercept
for(i in 2:p){
  X[,i] <- rnorm(n = N, mean = 0 , sd =1 ) 
}

# true parameters

beta <- runif(p,min=-1,max=1)
mu <- X %*% beta

sigma <- 1
#for repeated simulations
sigma <- runif(1, min = 1, max = 5)

alpha <- 1/sigma 
lambda <- exp(mu)/((log(2))^(1/alpha)) #attenzione alle diverse parametrizzazione della weibull
perc=0.2 #Percentuale dati censurati. Oppure 0.5 o 0.8
#perc=runif(1, min = 0.2, max = 0.8) #simulazioni ripetute
theta=(1/perc-1)^(1/alpha)*lambda


survt = rweibull(N, shape=alpha, scale = lambda) 
cent = rweibull(N, shape=alpha, scale = theta)

censured=survt>cent
delta <- as.logical(1-censured)

survt[delta==0] <- cent[delta==0] # censor survival time.



# count number of missing/censored survival times
n_miss <- N-sum(delta)

# data for censored subjects
y_m=survt[delta==0]
X_m=X[delta==0,]

# data for uncensored subjects
y_o=survt[delta==1]
X_o=X[delta==1,]
N_m = n_miss
N_o = N - n_miss

# initial values of parameters
beta0=rep(0.1, p) 
alpha0=0.1



###STAN

data_aft_NI <-list(p=p,
                   N_m = N_m,
                   X_m=as.matrix(X_m),
                   y_m=y_m,
                   N_o = N_o,
                   X_o=as.matrix(X_o),
                   y_o=y_o)


#initialization of the parameters

inits_NI <- function() 
{
  list(
    beta = beta0, 
    alpha=alpha0
  )
}

# run stan model
t<- system.time(
  {
    AFT_NI <- stan(file = "AFT_Non_Informative.stan", 
                   data = data_aft_NI,
                   chains = 1, 
                   iter = 10000, 
                   warmup = 5000, 
                   thin= 2, 
                   seed = 42, 
                   init = inits_NI,
                   algorithm = 'NUTS')
  })
