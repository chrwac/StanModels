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
  
  real conc0[num_traces];
  real conc0_oligos[num_traces];
  
  real<lower=0> k_act;
  real<lower=0> k_on;

  real<lower=0> bgs[num_traces];
  real<lower=0> amps[num_traces];
}

transformed data
{
  real x_r[0];
  int x_i[0];
  real thetas_dim[2];
  thetas_dim[1] = k_act;
  thetas_dim[2] = k_on;
}

model
{
}

generated quantities
{
  real intensities[num_times,num_traces];
  real conc_species[num_times,4];
  real initial_concs[4];
  initial_concs[3] = 0.0;
  initial_concs[4] = 0.0;
  for(i in 1:num_traces)
  {
    initial_concs[1] = conc0[i];
    initial_concs[2] = conc0_oligos[i];
    conc_species = integrate_ode_rk45(KineticSysDimWithActivation,initial_concs,0.0,times,thetas_dim,x_r,x_i);
    for(j in 1:num_times)
    {
      intensities[j,i] = bgs[i] + amps[i]*conc_species[j,4];
    }
  }
}
