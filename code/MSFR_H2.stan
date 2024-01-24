functions {
  vector MSFRm(real t, vector y, vector theta, int m) {
    vector[m] dydt;
    real denominator;
    denominator = 1.0;
    for(i in 1:m) denominator = denominator+(theta[i]*theta[m+i]*y[i]);
    for(i in 1:m) dydt[i] = -theta[i]*y[i] / denominator;
    for(i in 1:m){ if(y[i]<1e-6) dydt[i] = 0.0; }
    return dydt;
  }
}

data {
  int n;
  int m;
  array[n,m] int N_0;
  array[n,m] int N_e;
  array[n] real Time;
  array[n] int Replaced;
}

transformed data {
  matrix[n,m] N_0_real;
  int nLL=0;
  for(i in 1:n)
    for(j in 1:m)
      N_0_real[i,j] = 1.0*N_0[i,j];
  for(i in 1:n)
    for(j in 1:m)
      if(N_0[i,j]>0)
        nLL=nLL+1;
}

parameters{
  vector<lower=0>[m] a;
  vector<lower=0>[m] h;
}

model{
  array[1] vector[m] N_p; // predictions
  vector[2*m] pars;       // parameters for ODE
  
  // priors
  a ~ exponential(1);
  h ~ exponential(10);

  for (i in 1:m){
    pars[i] = a[i];
    pars[m+i] = h[i];
  }
 
  for(i in 1:n){
    if(Replaced[i]>0){
      N_p = { MSFRm(0.0, N_0_real[i, ]', pars, m) };
      for(j in 1:m)
        if(N_0[i,j]>0)   
          N_e[i,j] ~ poisson(-Time[i]*N_p[1][j]);
    } else {
      N_p = ode_rk45(MSFRm, N_0_real[i, ]', 0.0, {Time[i]}, pars, m);
      for(j in 1:m)
        if(N_0[i,j]>0)   
          N_e[i,j] ~ binomial(N_0[i,j], (N_0[i,j]-N_p[1][j])/N_0[i,j] );
    }
  }
}

generated quantities{
  array[nLL] real log_lik;{
    int k=0;
    array[1] vector[m] N_p; // predictions
    vector[2*m] pars;       // parameters for ODE
    for (i in 1:m){
      pars[i] = a[i];
      pars[m+i] = h[i];
    }
    for(i in 1:n){
      if(Replaced[i]>0){
        N_p = { MSFRm(0.0, N_0_real[i, ]', pars, m) };
        for(j in 1:m)
          if(N_0[i,j]>0){
            k=k+1;
            log_lik[k] = poisson_lpmf(N_e[i,j] | -Time[i]*N_p[1][j] );
          }
      } else {
        N_p = ode_rk45(MSFRm, N_0_real[i, ]', 0.0, {Time[i]}, pars, m);
        for(j in 1:m)
          if(N_0[i,j]>0){
            k=k+1;
            log_lik[k] = binomial_lpmf(N_e[i,j] | N_0[i,j], (N_0[i,j]-N_p[1][j])/N_0[i,j] );
          }
      }
    }
  }
}
