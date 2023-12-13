rm(list=ls())
library("rstan")

# Holling 2 --------------------------------------------------------------------

code.H2 = "
functions {
  real[] MSFRm(real t, real[] y, real[] theta, real[] x_r, int[] m) {
    real dydt[m[1]];
    real denominator;
    denominator = 1.0;
    for(i in 1:m[1]) denominator = denominator+(theta[i]*theta[m[1]+i]*y[i]);
    for(i in 1:m[1]) dydt[i] = -theta[i]*y[i] / denominator;
    for(i in 1:m[1]){ if(y[i]<1e-6) dydt[i] = 0.0; }
    return dydt;
  }
}

data {
  int n;
  int m;
  int N_0[n,m];
  int N_e[n,m];
  real Time[n];
  int Replaced[n];
}

transformed data {
  real x_r[0];
  real N_0_real[n,m];
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
  real N_p[1,m];     // predictions
  real N_p_rep[m];   // predictions
  real pars[2*m];    // parameters for ODE

  a ~ exponential(1);
  h ~ exponential(10);

  for (i in 1:m){
    pars[i] = a[i];
    pars[m+i] = h[i];
  }
 
  for(i in 1:n){
    if(Replaced[i]>0){
      N_p_rep = MSFRm(0.0,
                      N_0_real[i, ],
                      pars,
                      x_r,
                      {m} );
      for(j in 1:m)
        if(N_0[i,j]>0)   
          N_e[i,j] ~ poisson(-Time[i]*N_p_rep[j]);
    } else {
      N_p = integrate_ode_rk45(MSFRm, 
                               N_0_real[i, ],
                               0.0,
                               {Time[i]},
                               pars,
                               x_r, 
                               {m} );
      for(j in 1:m)
        if(N_0[i,j]>0)   
          N_e[i,j] ~ binomial(N_0[i,j], (N_0[i,j]-N_p[1,j])/N_0[i,j] );
    }
  }
}

generated quantities{
  real log_lik[nLL];{
    int k=0;
    real N_p[1,m];     // predictions
    real N_p_rep[m];   // predictions
    real pars[2*m];    // parameters for ODE
    for (i in 1:m){
      pars[i] = a[i];
      pars[m+i] = h[i];
    }
    for(i in 1:n){
      if(Replaced[i]>0){
        N_p_rep = MSFRm(0.0,
                        N_0_real[i, ],
                        pars,
                        x_r,
                        {m} );
        for(j in 1:m)
          if(N_0[i,j]>0){
            k=k+1;
            log_lik[k] = poisson_lpmf(N_e[i,j] | -Time[i]*N_p_rep[j] );
          }
      } else {
        N_p = integrate_ode_rk45(MSFRm, 
                                 N_0_real[i, ],
                                 0.0,
                                 {Time[i]},
                                 pars,
                                 x_r, 
                                 {m} );
        for(j in 1:m)
          if(N_0[i,j]>0){
            k=k+1;
            log_lik[k] = binomial_lpmf(N_e[i,j] | N_0[i,j], (N_0[i,j]-N_p[1,j])/N_0[i,j] );
          }
      }
    }
  }
}
"

model.H2 = stan_model(model_code = code.H2)
save(model.H2, file="stan_model_H2_mixed.RData")
remove(model.H2)

print("Holling 2 OK")


# Holling 3 --------------------------------------------------------------------

code.H3 = "
functions {
  real[] MSFRm(real t, real[] y, real[] theta, real[] x_r, int[] m) {
    real dydt[m[1]];
    real denominator;
    denominator = 1.0;
    for(i in 1:m[1]) denominator = denominator+(theta[i]*theta[m[1]+i]*y[i]^(1.0+theta[1+2*m[1]]));
    for(i in 1:m[1]) dydt[i] = -theta[i]*y[i]^(1.0+theta[1+2*m[1]]) / denominator;
    for(i in 1:m[1]){ if(y[i]<1e-6) dydt[i] = 0.0; }
    return dydt;
  }
}

data {
  int n;
  int m;
  int N_0[n,m];
  int N_e[n,m];
  real Time[n];
  int Replaced[n];
}

transformed data {
  real x_r[0];
  real N_0_real[n,m];
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
  real<lower=0> q;
}

model{
  real N_p[1,m];     // predictions
  real N_p_rep[m];   // predictions
  real pars[1+2*m];  // parameters for ODE

  a ~ exponential(1);
  h ~ exponential(10);
  q ~ exponential(1);
 
  for (i in 1:m){
    pars[i] = a[i];
    pars[m+i] = h[i];
  }
  pars[1+2*m] = q;
 
  for(i in 1:n){
    if(Replaced[i]>0){
      N_p_rep = MSFRm(0.0,
                      N_0_real[i, ],
                      pars,
                      x_r,
                      {m} );
      for(j in 1:m)
        if(N_0[i,j]>0)   
          N_e[i,j] ~ poisson(-Time[i]*N_p_rep[j]);
    } else {
      N_p = integrate_ode_rk45(MSFRm, 
                               N_0_real[i, ],
                               0.0,
                               {Time[i]},
                               pars,
                               x_r, 
                               {m} );
      for(j in 1:m)
        if(N_0[i,j]>0)   
          N_e[i,j] ~ binomial(N_0[i,j], (N_0[i,j]-N_p[1,j])/N_0[i,j] );
    }
  }
}

generated quantities{
  real log_lik[nLL];{
    int k=0;
    real N_p[1,m];     // predictions
    real N_p_rep[m];   // predictions
    real pars[1+2*m];  // parameters for ODE
    for (i in 1:m){
      pars[i] = a[i];
      pars[m+i] = h[i];
    }
    pars[1+2*m] = q;
    for(i in 1:n){
      if(Replaced[i]>0){
        N_p_rep = MSFRm(0.0,
                        N_0_real[i, ],
                        pars,
                        x_r,
                        {m} );
        for(j in 1:m)
          if(N_0[i,j]>0){
            k=k+1;
            log_lik[k] = poisson_lpmf(N_e[i,j] | -Time[i]*N_p_rep[j] );
          }
      } else {
        N_p = integrate_ode_rk45(MSFRm, 
                                 N_0_real[i, ],
                                 0.0,
                                 {Time[i]},
                                 pars,
                                 x_r, 
                                 {m} );
        for(j in 1:m)
          if(N_0[i,j]>0){
            k=k+1;
            log_lik[k] = binomial_lpmf(N_e[i,j] | N_0[i,j], (N_0[i,j]-N_p[1,j])/N_0[i,j] );
          }
      }
    }
  }
}
"

model.H3 = stan_model(model_code = code.H3)
save(model.H3, file="stan_model_H3_mixed.RData")
remove(model.H3)

print("Holling 3 OK")


# Yodzis  ---------------------------------------------------------------------

code.Yo = "
functions {
  real[] MSFRm(real t, real[] y, real[] theta, real[] x_r, int[] m) {
    real dydt[m[1]];
    real numerator[m[1]];
    real denominator = 0.0;
    for(i in 1:m[1]) numerator[i] = theta[i]*y[i]^theta[3*m[1]+2];
    for(i in 1:m[1]) denominator = denominator + theta[2*m[1]+i]*y[i]^theta[3*m[1]+1] + theta[m[1]+i]*numerator[i];
    for(i in 1:m[1]) dydt[i] = - numerator[i] / denominator;
    for(i in 1:m[1]){ if(y[i]<1e-6) dydt[i] = 0.0; }
    return dydt;
  }
}

data {
  int n;
  int m;
  int N_0[n,m];
  int N_e[n,m];
  real Time[n];
  int Replaced[n];
}

transformed data {
  real x_r[0];
  real N_0_real[n,m];
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
  simplex[m] w;
  real<lower=0.0> r;
}

model{
  real N_p[1,m];     // predictions
  real N_p_rep[m];   // predictions
  real pars[3*m+2];  // parameters for ODE
  vector[m] w_prior; // m-dim parametrisation for dirichlet

  for (i in 1:m){
    pars[i] = a[i]*w[i];
    pars[m+i] = h[i];
    pars[2*m+i] = w[i];
  }
  pars[3*m+1] = r;
  pars[3*m+2] = 1.0+r; // directly hand over sum

  for (i in 1:m){
    w_prior[i] = 2.0;
  }

  a ~ exponential(1);
  h ~ exponential(10);
  w ~ dirichlet(w_prior);
  r ~ exponential(1);

  for(i in 1:n){
    if(Replaced[i]>0){
      N_p_rep = MSFRm(0.0,
                      N_0_real[i, ],
                      pars,
                      x_r,
                      {m} );
      for(j in 1:m)
        if(N_0[i,j]>0)   
          N_e[i,j] ~ poisson(-Time[i]*N_p_rep[j]);
    } else {
      N_p = integrate_ode_rk45(MSFRm, 
                               N_0_real[i, ],
                               0.0,
                               {Time[i]},
                               pars,
                               x_r, 
                               {m} );
      for(j in 1:m)
        if(N_0[i,j]>0)   
          N_e[i,j] ~ binomial(N_0[i,j], (N_0[i,j]-N_p[1,j])/N_0[i,j] );
    }
  }
}

generated quantities{
  real log_lik[nLL];{
    int k=0;
    real N_p[1,m];     // predictions
    real N_p_rep[m];   // predictions
    real pars[3*m+2];  // parameters for ODE
    for (i in 1:m){
      pars[i] = a[i]*w[i];
      pars[m+i] = h[i];
      pars[2*m+i] = w[i];
    }
    pars[3*m+1] = r;
    pars[3*m+2] = 1.0+r; // directly hand over sum
    for(i in 1:n){
      if(Replaced[i]>0){
        N_p_rep = MSFRm(0.0,
                        N_0_real[i, ],
                        pars,
                        x_r,
                        {m} );
        for(j in 1:m)
          if(N_0[i,j]>0){
            k=k+1;
            log_lik[k] = poisson_lpmf(N_e[i,j] | -Time[i]*N_p_rep[j] );
          }
      } else {
        N_p = integrate_ode_rk45(MSFRm, 
                                 N_0_real[i, ],
                                 0.0,
                                 {Time[i]},
                                 pars,
                                 x_r, 
                                 {m} );
        for(j in 1:m)
          if(N_0[i,j]>0){
            k=k+1;
            log_lik[k] = binomial_lpmf(N_e[i,j] | N_0[i,j], (N_0[i,j]-N_p[1,j])/N_0[i,j] );
          }
      }
    }
  }
}
"

model.Yo = stan_model(model_code = code.Yo)
save(model.Yo, file="stan_model_Yo_mixed.RData")
remove(model.Yo)

print("Yodzis OK")


# Generalized switching --------------------------------------------------------

code.Gen = "
functions {
  real[] MSFRm(real t, real[] y, real[] theta, real[] x_r, int[] m) {
    real dydt[m[1]];
    real numerator[m[1]];
    real denominator = 0.0;
    for(i in 1:m[1]) numerator[i] = theta[i]*y[i]^theta[3*m[1]+2];
    for(i in 1:m[1]) denominator = denominator + theta[2*m[1]+i]*y[i]^theta[3*m[1]+1] + theta[m[1]+i]*numerator[i];
    for(i in 1:m[1]) dydt[i] = - numerator[i] / denominator;
    for(i in 1:m[1]){ if(y[i]<1e-6) dydt[i] = 0.0; }
    return dydt;
  }
}

data {
  int n;
  int m;
  int N_0[n,m];
  int N_e[n,m];
  real Time[n];
  int Replaced[n];
}

transformed data {
  real x_r[0];
  real N_0_real[n,m];
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
  simplex[m] w;
  real<lower=0.0> q;
  real<lower=0.0> r;
}

model{
  real N_p[1,m];     // predictions
  real N_p_rep[m];   // predictions
  real pars[3*m+2];  // parameters for ODE
  vector[m] w_prior; // m-dim parametrisation for dirichlet

  for (i in 1:m){
    pars[i] = a[i]*w[i];
    pars[m+i] = h[i];
    pars[2*m+i] = w[i];
  }
  pars[3*m+1] = r;       // exponent for Yodzis part
  pars[3*m+2] = 1.0+q+r; // directly hand over sum

  for (i in 1:m){
    w_prior[i] = 2.0;
  }
  a ~ exponential(0.1);
  h ~ exponential(1);
  w ~ dirichlet(w_prior);
  q ~ exponential(1);
  r ~ exponential(1);

  for(i in 1:n){
    if(Replaced[i]>0){
      N_p_rep = MSFRm(0.0,
                      N_0_real[i, ],
                      pars,
                      x_r,
                      {m} );
      for(j in 1:m)
        if(N_0[i,j]>0)   
          N_e[i,j] ~ poisson(-Time[i]*N_p_rep[j]);
    } else {
      N_p = integrate_ode_rk45(MSFRm, 
                               N_0_real[i, ],
                               0.0,
                               {Time[i]},
                               pars,
                               x_r, 
                               {m} );
      for(j in 1:m)
        if(N_0[i,j]>0)   
          N_e[i,j] ~ binomial(N_0[i,j], (N_0[i,j]-N_p[1,j])/N_0[i,j] );
    }
  }
}

generated quantities{
  real log_lik[nLL];{
    int k=0;
    real N_p[1,m];     // predictions
    real N_p_rep[m];   // predictions
    real pars[3*m+2];  // parameters for ODE
    for (i in 1:m){
      pars[i] = a[i]*w[i];
      pars[m+i] = h[i];
      pars[2*m+i] = w[i];
    }
    pars[3*m+1] = r;       // exponent for Yodzis part
    pars[3*m+2] = 1.0+q+r; // directly hand over sum
    for(i in 1:n){
      if(Replaced[i]>0){
        N_p_rep = MSFRm(0.0,
                        N_0_real[i, ],
                        pars,
                        x_r,
                        {m} );
        for(j in 1:m)
          if(N_0[i,j]>0){
            k=k+1;
            log_lik[k] = poisson_lpmf(N_e[i,j] | -Time[i]*N_p_rep[j] );
          }
      } else {
        N_p = integrate_ode_rk45(MSFRm, 
                                 N_0_real[i, ],
                                 0.0,
                                 {Time[i]},
                                 pars,
                                 x_r, 
                                 {m} );
        for(j in 1:m)
          if(N_0[i,j]>0){
            k=k+1;
            log_lik[k] = binomial_lpmf(N_e[i,j] | N_0[i,j], (N_0[i,j]-N_p[1,j])/N_0[i,j] );
          }
      }
    }
  }
}
"

model.Gen = stan_model(model_code = code.Gen)
save(model.Gen, file="stan_model_Gen_mixed.RData")
remove(model.Gen)

print("New OK")


