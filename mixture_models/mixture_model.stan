//#################################################
//###  MIXTURE MODEL ###########
//#################################################




data {
  int <lower=1> N;
  vector[N] Y;
  int H;
}
parameters {
  ordered[H] mu;
  vector<lower = 0>[H] sigma2;
  simplex[H] theta;
real <lower=0> v;
}

model {
  vector[H] contributions;
  // priors
  mu ~ normal(0, sqrt(v));
  sigma2 ~ inv_gamma(1, 1);
  theta ~ dirichlet(rep_vector(1.0, H));
  v ~ inv_gamma(1, 1);
  
  // likelihood
  for(i in 1:N) {
    for(k in 1:H) {
      contributions[k] = log(theta[k]) + normal_lpdf(Y[i] | mu[k], pow(sigma2[k], 0.5));
    }
    target += log_sum_exp(contributions);
  }
}