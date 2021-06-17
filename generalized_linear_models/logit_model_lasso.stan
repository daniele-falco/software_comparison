//#################################################
//###  PROBIT REGRESSION ###########
//#################################################

///////////////////////// DATA /////////////////////////////////////
data {
	int<lower = 0> N;       // number of data
	int<lower = 0> p;       // number of covariates (without the intercept)
	int<lower = 0, upper = 1> Y[N];   // response vector
	matrix[N, p] X;   		  // design matrix

}

//////////////////// PARAMETERS /////////////////////////////////
parameters {
	vector[p] beta;        // regression coefficients 
	 real <lower=0> lambda2; 	        
	
}


////////////////// MODEL ////////////////////////
model {

// Likelihood     
	Y ~ bernoulli_logit( X*beta );
// Prior
	 

	for (j in 1:p) 
	{
	 	beta[j] ~ double_exponential(0.0, 1.0 / (pow(lambda2,0.5)));
	}
lambda2 ~ exponential(0.1);
	

}


