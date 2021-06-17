
//

data {
  int<lower=0> N; // number of data
  vector[N] y;    // response
  int<lower=0> p; //number of regressors
  matrix[N, p] X; //matrix of covariates
  

}


parameters {
  vector[p] beta;       //regressors paramters (column vector)

  real<lower=0> sigma;
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
  y ~ normal(mu, sigma);

for (j in 1:p){
	beta[j]~ normal(0,100);}
	sigma~cauchy(0, 2.5)T[0,];


}
