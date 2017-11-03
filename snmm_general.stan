data
{
  int<lower=1> N; // number of data points
  int<lower=1> K; // number of mixtures/components
  real y[N]; // data points
  vector<lower=0>[K] alphas; // priors for the mixture components
  real<lower=0> mus0[K]; // priors for the mixture means
  real<lower=0> std0; // the hyperparameter: standard-deviation for the means...
}

parameters
{
  simplex[K] pis; 
  real locations[K];
  real<lower=0> scales[K];
  real shapes[K];
}

model
{
  real ps[K];
  locations ~ normal(mus0,std0);
  scales ~ lognormal(-3,2);
  shapes ~ normal(0,0.5);
  pis ~ dirichlet(alphas);

  for(n in 1:N)
  {
    for(i in 1:K)
    {
      ps[i] <- log(pis[i]) + skew_normal_log(y[n],locations[i],scales[i],shapes[i]);
    }
    increment_log_prob(log_sum_exp(ps));
  }
}