This repository contains the code use for the comparison of three Bayesian MCMC software: JAGS, Stan and NIMBLE.
Each folder is related to a specific class of models.
We considered linear models, generalized linear models (i.e. logistic regression), mixture models and accelerated failure time models.
Each folder contains a .R file for each considered prior distribution.
The .R file contains the generation of the data, the model specification (unless those related to Stan which are contained in a separate .STAN file), and the generation of the MCMC chains.
