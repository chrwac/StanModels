functions
{
  real[] KineticSysDimWithActivation(real t,real[] y,real[] theta,real[] x_r,int[] x_i)
  {
    real dydt[4];
    dydt[1] = -theta[1]*y[1]*y[2];
    dydt[2] = -theta[1]*y[1]*y[2];
    dydt[3] = theta[1]*y[1]*y[2] - theta[2]*y[3]*y[3];
    dydt[4] = theta[2]*y[3]*y[3];
    return dydt;
  }
}

data
{
  int<lower=1> num_traces; // how many different traces (conc0,conc0_oligos) - pairs ?
  int<lower=1> num_times; //
  
  real<lower=0> times[num_times];
  real intensities[num_traces,num_times];
  
  real conc0[num_traces];
  real conc0_oligos[num_traces];
}

transformed data
{
  real x_r[0];
  int x_i[0];
  
}

parameters
{
  real<lower=0> k_act;
  real<lower=0> k_on;
  real<lower=0> sigma;
  real<lower=0> bgs[num_traces];
  real<lower=0> amps[num_traces];
}

model
{
  real conc_species[num_times,4];
  real intens_nf[num_traces,num_times];
  
  real concs0[4];
  real thetas_dim[2];
  //real times_with_offest[num_times];
  
  k_act ~ lognormal(-8,2);
  k_on ~ lognormal(-8,2);
  sigma ~ lognormal(3,2);
  
  bgs ~ normal(1000,200);
  amps ~ normal(780,10);
  
  thetas_dim[1] = k_act;
  thetas_dim[2] = k_on;
  
  concs0[3] = 0;
  concs0[4] = 0;
  
  for(i in 1:num_traces)
  {
    concs0[1] =  conc0[i];
    concs0[2] = conc0_oligos[i];
    conc_species = integrate_ode_rk45(KineticSysDimWithActivation,concs0,0.0,times,thetas_dim,x_r,x_i);
    for(j in 1:num_times)
    {
      intens_nf[i,j] = bgs[i] + amps[i]*conc_species[j,4];
      intensities[i,j] ~ normal(intens_nf[i,j],sigma);
    }
  }
}
