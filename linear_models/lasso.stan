
//

data {
  int<lower=0> N; // number of data
  vector[N] y;    // response
  int<lower=0> p; //number of regressors
  matrix[N, p] X; //matrix of covariates
  real sigma0;
  real nu0;


}


parameters {
  vector[p] beta;       //regressors paramters (column vector)
real lambda2;
  real<lower=0> sigma2;
}

transformed parameters 
{
	vector[N] mu;         // mean 

	for(i in 1:N){
    mu[i] =  X[i,] * beta;
	}



}

model {
  //likelihood:
  y ~ normal(mu, pow(sigma2, 0.5));
for (j in 1:p) 
	{
	 	beta[j] ~ double_exponential(0.0, 1.0 / (pow(lambda2,0.5)));
	}
	sigma2 ~ inv_gamma(nu0/2,nu0*sigma0*sigma0/2);
	lambda2 ~ exponential(0.1);


}
