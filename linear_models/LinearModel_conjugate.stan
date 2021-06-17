
//

data {
  int<lower=0> N; // number of data
  vector[N] y;    // response
  int<lower=0> p; //number of regressors
  matrix[N, p] X; //matrix of covariates
  
  
  matrix[p, p] B0;
  vector[p] b0; 
   
  real sigma0;
	real nu0; 
}


parameters {
  vector[p] beta;       //regressors paramters (column vector)
	real sigma2;
  
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
  
  //Prior:
  beta ~ multi_normal(b0,  sigma2*B0);
	sigma2 ~ inv_gamma(nu0/2,nu0*sigma0*sigma0/2 );

}

