rm(list=ls())
setwd("~/Nextcloud/MSFR/Tutorial")

library("rstan")
library("loo")
library("coda")
library("deSolve")
library("BayesianTools")

set.seed(100)

fit.new = TRUE

# compile stan models ----------------------------------------------------------
# source("compilation_loo_mixed.R")

# source functions -------------------------------------------------------------
# for generating design and ODE functions for desolve 
source("utils_datasim.R")
source("desolve_models.R")

# simulate dataset -------------------------------------------------------------

m = 2 # number of species
n = 256 # number of feeding trials: 64, 128, 256, 512, 1024
Nmax = 200 # maximum number of offered prey per species

strait = 0.1 # between [0.0, 0.4] in the manuscript

DOE = doe.gridlog # choices? 3d data?
MSFR = Yodzis

df.all = DOE(n, Nmax) |> as.data.frame()

colnames(df.all) = c("Prey1", "Prey2", "Prey1.Eaten", "Prey2.Eaten" )
df.all$Replaced = TRUE
df.all$Time = 1.0
head(df.all)

# mean traits
Fmax = runif(m, 50, 100)
Nhalf = runif(m, 25, 75)
w = c(0.5,0.5)
r = 1.0

# pars.true = list(a=Fmax/Nhalf, h=1/Fmax) # Holling2
# pars.true = list(a=Fmax/Nhalf^(1+q), h=1/Fmax, q=q) # Holling3
pars.true = list(a=Fmax/Nhalf, h=1/Fmax, w=w, r=r) # Yodzis
# pars.true = list(a=Fmax/Nhalf^(1+q), h=1/Fmax, w=w, q=q, r=r) # New

pars.true = append(pars.true, list(m=m, n=n, strait=strait, DOE=DOE))
pars.individual = pars.true
ai = rep(0,m)
hi = rep(0,m)

for(i in 1:n){ # n feeding trials 
  # traits
  for(j in 1:m){ # each species
    ai[j] = max(rnorm(1, pars.true$a[j], strait*pars.true$a[j]), 
                pars.true$a[j]*(1-2*strait))
    hi[j] = max(rnorm(1, pars.true$h[j], strait*pars.true$h[j]), 
                pars.true$h[j]*(1-2*strait))
  }
  pars.individual$a = ai
  pars.individual$h = hi
  if(df.all$Replaced[i]==0){
    # simulation
    out = ode(y = as.numeric(df.all[i, 1:m]),
              times = c(0,df.all$Time[i]),
              func = MSFR,
              parms = pars.individual)
    # integer eaten prey
    for(j in 1:m){
      df.all[i,m+j] = df.all[i,j]-out[2,j+1] # eaten prey
      if(df.all[i,m+j]>0){ # draw from binomial
        df.all[i,m+j] = rbinom(1, size=df.all[i,j], prob=df.all[i,m+j]/df.all[i,j])
      }
    }
  } else {
    # simulation
    out = (-1)*unlist(MSFR(N=as.numeric(df.all[i, 1:m]),
                           parms = pars.individual))
    # integer eaten prey
    for(j in 1:m){
      df.all[i,m+j] = df.all$Time[i] * out[j] # eaten prey
      if(df.all[i,m+j]>0){ # draw from poisson
        df.all[i,m+j] = rpois(1, df.all[i,m+j])
      }
    }
  }
}

head(df.all)
df.single.1 = subset(df.all, Prey2==0)
df.single.2 = subset(df.all, Prey1==0)
df.multi = subset(df.all, Prey1>0 & Prey2>0)

par(mfrow=c(2,2))
plot(df.all[, 1:2], pch=16, xlab="Prey 1", ylab="Prey 2", main="Experimental design")
plot(df.single.1$Prey1, df.single.1$Prey1.Eaten, 
     xlim=c(0,Nmax), ylim=c(0,max(df.all[, 3:4])),
     xlab="Prey 1 offered", ylab="Prey 1 eaten", main="Single-species feeding: Prey 1")
plot(df.single.2$Prey2, df.single.2$Prey2.Eaten, 
     xlim=c(0,Nmax), ylim=c(0,max(df.all[, 3:4])),
     xlab="Prey 2 offered", ylab="Prey 2 eaten", main="Single-species feeding: Prey 2")
plot(df.multi$Prey1/(df.multi$Prey1+df.multi$Prey2),
     df.multi$Prey1.Eaten/(df.multi$Prey1.Eaten+df.multi$Prey2.Eaten),
     xlim=c(0,1), ylim=c(0,1),
     xlab="Proportion prey 1 offered", ylab="Proportion prey 1 eaten", main="Multi-species preference")
abline(0,1)

data.stan = list(n=pars.true$n,
                 m=pars.true$m,
                 N_0=df.all[, 1:pars.true$m], 
                 N_e=df.all[, (pars.true$m+1):(2*pars.true$m)],
                 Replaced = df.all$Replaced,
                 Time = df.all$Time
)

# fit H2 -----------------------------------------------------------------------

load("stan_model_H2_mixed.RData")

chains = 3

# initial values for sampling
init = rep(list(list(a = rep(1,pars.true$m),
                     h = rep(0.01,pars.true$m)
)),chains)

if(fit.new){
  fit1 = sampling(
    model.H2,
    data = data.stan,
    chains = chains,
    warmup = 1000,
    iter = 3000,
    cores = 3,
    init = init,
    refresh = 100
  )
  save(fit1, file=("stan_fit_H2.RData"))
} else{
  load(file=("stan_fit_H2.RData"))
}

print(fit1, pars=c("a","h"), digits=3, probs=c(0.05, 0.5, 0.95))

# fit H3 -----------------------------------------------------------------------

load("stan_model_H3_mixed.RData")

chains = 3

# initial values for sampling
init = rep(list(list(a = rep(0.1,pars.true$m),
                     h = rep(0.01,pars.true$m),
                     q = 1.0
)),chains)

if(fit.new){
  fit2 = sampling(
    model.H3,
    data = data.stan,
    chains = chains,
    warmup = 1000,
    iter = 3000,
    cores = 3,
    init = init,
    refresh = 100
  )
  save(fit2, file=("stan_fit_H3.RData"))
} else{
  load(file=("stan_fit_H3.RData"))
}

print(fit2, pars=c("a","h","q"), digits=3, probs=c(0.05, 0.5, 0.95))

# fit Yo -----------------------------------------------------------------------

load("stan_model_Yo_mixed.RData")

chains = 3

# initial values for sampling
init = rep(list(list(a = rep(1,pars.true$m),
                     h = rep(0.01,pars.true$m),
                     r = 1.0
)),chains)

if(fit.new){
  fit3 = sampling(
    model.Yo,
    data = data.stan,
    chains = chains,
    warmup = 1000,
    iter = 3000,
    cores = 3,
    init = init,
    refresh = 100
  )
  save(fit3, file=("stan_fit_Yo.RData"))
} else{
  load(file=("stan_fit_Yo.RData"))
}

print(fit3, pars=c("a","h","w","r"), digits=3, probs=c(0.05, 0.5, 0.95))

# fit New ----------------------------------------------------------------------

load("stan_model_New_mixed.RData")

chains = 3

# initial values for sampling
init = rep(list(list(a = rep(0.1,pars.true$m),
                     h = rep(0.01,pars.true$m),
                     q = 1.0,
                     r = 1.0
)),chains)

if(fit.new){
  fit4 = sampling(
    model.New,
    data = data.stan,
    chains = chains,
    warmup = 1000,
    iter = 3000,
    cores = 3,
    init = init,
    refresh = 100
  )
  save(fit4, file=("stan_fit_New.RData"))
} else{
  load(file=("stan_fit_New.RData"))
}

print(fit4, pars=c("a","h","w","q","r"), digits=3, probs=c(0.05, 0.5, 0.95))

# model comparison -------------------------------------------------------------

lik1 = extract_log_lik(fit1, merge_chains=FALSE)
reff1 = relative_eff(exp(lik1), cores=4) 
loo1 = loo(lik1, r_eff=reff1, cores=4)

lik2 = extract_log_lik(fit2, merge_chains=FALSE)
reff2 = relative_eff(exp(lik2), cores=4) 
loo2 = loo(lik2, r_eff=reff2, cores=4)

lik3 = extract_log_lik(fit3, merge_chains=FALSE)
reff3 = relative_eff(exp(lik3), cores=4) 
loo3 = loo(lik3, r_eff=reff3, cores=4)

lik4 = extract_log_lik(fit4, merge_chains=FALSE)
reff4 = relative_eff(exp(lik4), cores=4) 
loo4 = loo(lik4, r_eff=reff4, cores=4)

loo_compare(loo1, loo2, loo3, loo4)

# posterior predictive checks --------------------------------------------------

MSFR = Yodzis

n.post = 200
n.pred = 20

post = as.matrix(fit3)[, 1:7]
set.seed(100)
post = post[sample(1:nrow(post),n.post), ]
head(post)

par(mfrow=c(1,1))

# prediction single 1 ----------------------------------------------------------

N10.pred = seq(0.1, Nmax, length.out=n.pred)
N1E.pred = matrix(NA, nrow=n.post, ncol=n.pred)

replaced = TRUE
time = 1.0

for(k in 1:n.post){
  params = list(a=post[k, 1:2], 
                h=post[k, 3:4],
                w=post[k, 5:6],
                r=post[k, 7],
                m=2)
  for(i in 1:n.pred){
    if(replaced==FALSE){
      out = ode(y = c(N=c(N10.pred[i], 0.0)), times = c(0,time), func = MSFR, parms = params) 
      N1E.pred[k,i] = N10.pred[i]-out[2,2]
    } else{
      out = (-1)*unlist(MSFR(N=c(N10.pred[i], 0.0), parms=params))
      N1E.pred[k,i] = out[1]*time
    }
  }
}

N1E.pred.qs = apply(N1E.pred, 2, function(x) quantile(x, probs=c(0.05, 0.5, 0.95)))

plot(df.single.1$Prey1, df.single.1$Prey1.Eaten, 
     xlim=c(0,Nmax), ylim=c(0,max(df.all[, 3:4])),
     xlab="Prey 1 offered", ylab="Prey 1 eaten", main="Single-species feeding: Prey 1")
polygon( c(N10.pred, rev(N10.pred)), 
         c(N1E.pred.qs[1, ], rev(N1E.pred.qs[3, ])), 
         border=NA,
         col=adjustcolor("blue", alpha.f=0.2)
)
lines(N10.pred, N1E.pred.qs[2, ], col="blue", lwd=2)

# prediction single 2 ----------------------------------------------------------

N20.pred = seq(0.1, Nmax, length.out=n.pred)
N2E.pred = matrix(NA, nrow=n.post, ncol=n.pred)

replaced = TRUE
time = 1.0

for(k in 1:n.post){
  params = list(a=post[k, 1:2], 
                h=post[k, 3:4],
                w=post[k, 5:6],
                r=post[k, 7],
                m=2)
  for(i in 1:n.pred){
    if(replaced==FALSE){
      out = ode(y = c(N=c(0.0, N20.pred[i])), times = c(0,time), func = MSFR, parms = params) 
      N2E.pred[k,i] = N20.pred[i]-out[2,3]
    } else{
      out = (-1)*unlist(MSFR(N=c(0.0, N20.pred[i]), parms=params))
      N2E.pred[k,i] = out[2]*time
    }
  }
}

N2E.pred.qs = apply(N2E.pred, 2, function(x) quantile(x, probs=c(0.05, 0.5, 0.95)))

plot(df.single.2$Prey2, df.single.2$Prey2.Eaten, 
     xlim=c(0,Nmax), ylim=c(0,max(df.all[, 3:4])),
     xlab="Prey 2 offered", ylab="Prey 2 eaten", main="Single-species feeding: Prey 2")
polygon( c(N20.pred, rev(N20.pred)), 
         c(N2E.pred.qs[1, ], rev(N2E.pred.qs[3, ])), 
         border=NA,
         col=adjustcolor("red", alpha.f=0.2)
)
lines(N20.pred, N2E.pred.qs[2, ], col="red", lwd=2)

# prediction mix ---------------------------------------------------------------

N0.total = 200

N10.pred = seq(0, N0.total, length.out=n.pred)
N20.pred = N0.total-N10.pred
N1E.pred = matrix(NA, nrow=n.post, ncol=n.pred)
N2E.pred = matrix(NA, nrow=n.post, ncol=n.pred)

replaced = TRUE
time = 1.0

for(k in 1:n.post){
  params = list(a=post[k, 1:2], 
                h=post[k, 3:4],
                w=post[k, 5:6],
                r=post[k, 7],
                m=2)
  for(i in 1:n.pred){
    if(replaced==FALSE){
      out = ode(y = c(N=c(N10.pred[i], N20.pred[i])), times = c(0,time), func = MSFR, parms = params) 
      N1E.pred[k,i] = N10.pred[i]-out[2,2]
      N2E.pred[k,i] = N20.pred[i]-out[2,3]
    }
    else{
      out = (-1)*unlist(MSFR(N=c(N10.pred[i], N20.pred[i]), parms=params))
      N1E.pred[k,i] = out[1]
      N2E.pred[k,i] = out[2]
    }
  }
}

# make relative 
N1E.pred.rel = matrix(NA, nrow=n.post, ncol=n.pred)
N2E.pred.rel = matrix(NA, nrow=n.post, ncol=n.pred)
for(k in 1:n.post){
  N1E.pred.rel[k, ] = N1E.pred[k, ]/(N1E.pred[k, ]+N2E.pred[k, ])
  N2E.pred.rel[k, ] = N2E.pred[k, ]/(N1E.pred[k, ]+N2E.pred[k, ])
}

N1E.pred.qs = apply(N1E.pred, 2, function(x) quantile(x, probs=c(0.05, 0.5, 0.95)))
N2E.pred.qs = apply(N2E.pred, 2, function(x) quantile(x, probs=c(0.05, 0.5, 0.95)))
N1E.pred.rel.qs = apply(N1E.pred.rel, 2, function(x) quantile(x, probs=c(0.05, 0.5, 0.95)))
N2E.pred.rel.qs = apply(N2E.pred.rel, 2, function(x) quantile(x, probs=c(0.05, 0.5, 0.95)))

plot(df.multi$Prey1/(df.multi$Prey1+df.multi$Prey2),
     df.multi$Prey1.Eaten/(df.multi$Prey1.Eaten+df.multi$Prey2.Eaten),
     xlim=c(0,1), ylim=c(0,1),
     xlab="Proportion prey 1 offered", ylab="Proportion prey 1 eaten", main="Multi-species preference")
abline(0,1)
polygon( c(N10.pred/N0.total, rev(N10.pred/N0.total)), 
         c(N1E.pred.rel.qs[1, ], rev(N1E.pred.rel.qs[3, ])), 
         border=NA,
         col=adjustcolor("blue", alpha.f=0.2)
)
lines(N10.pred/N0.total, N1E.pred.rel.qs[2, ], col="blue", lwd=2)

# obs vs pred mix ---------------------------------------------------------------

N1E.pred = matrix(NA, nrow=n.post, ncol=nrow(df.multi))
N2E.pred = matrix(NA, nrow=n.post, ncol=nrow(df.multi))

for(k in 1:n.post){
  params = list(a=post[k, 1:2], 
                h=post[k, 3:4],
                w=post[k, 5:6],
                r=post[k, 7],
                m=2)
  for(i in 1:nrow(df.multi)){
    if(df.all$Replaced[i]==FALSE){
      out = ode(y = c(N=c(df.multi$Prey1[i],df.multi$Prey2[i])), times = c(0,df.multi$Time[i]), func = MSFR, parms = params) 
      N1E.pred[k,i] = df.multi$Prey1[i]-out[2,2]
      N2E.pred[k,i] = df.multi$Prey2[i]-out[2,3]
    } else {
      out = (-1)*unlist(MSFR(N=c(df.multi$Prey1[i],df.multi$Prey2[i]), parms=params))
      N1E.pred[k,i] = out[1]*df.multi$Time[i]
      N2E.pred[k,i] = out[2]*df.multi$Time[i]
    }
  }
}

N1E.pred.qs = apply(N1E.pred, 2, function(x) quantile(x, probs=c(0.05, 0.5, 0.95)))
N2E.pred.qs = apply(N2E.pred, 2, function(x) quantile(x, probs=c(0.05, 0.5, 0.95)))

plot(NULL, 
     xlim=c(0,max(df.all[, 3:4])), ylim=c(0,max(df.all[, 3:4])),
     xlab="Predicted", ylab="Observed", main="Multi-species obs vs. pred")
points(N1E.pred.qs[2, ], df.multi$Prey1.Eaten, col="blue")
abline(0,1)
points(N2E.pred.qs[2, ], df.multi$Prey2.Eaten, col="red")
abline(0,1)

# Stan code --------------------------------------------------------------------
