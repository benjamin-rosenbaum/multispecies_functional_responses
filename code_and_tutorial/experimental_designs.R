
library("randtoolbox") # for Halton design

doe.random = function(n,Nmax){
  data = matrix(NA, nrow=n, ncol=4) 
  for(i in 1:2){
    data[, i] = sample(0:Nmax, size=n, replace=TRUE)
  }
  return(data)
}

doe.halton = function(n,Nmax){
  data = matrix(NA, nrow=n, ncol=4)
  dummy = halton(sample(1:1e6, size=1), dim=2) # halton of random size for init
  data[, 1:2] = round(Nmax*halton(n, dim=2, init=FALSE))
  return(data)
}

doe.triangle = function(n,Nmax){
  data = matrix(NA, nrow=n, ncol=4)
  n.triangle = floor(n/4) # divide 1/4 per single species, 1/2 for mixed
  triangle.axis = (1:n.triangle)*Nmax/n.triangle
  N0 = as.matrix( rbind( expand.grid(x=triangle.axis, y=0),
                         expand.grid(y=triangle.axis, x=0)))
  N1.0 = seq(from=0, to=Nmax, length.out=2*n.triangle)
  N2.0 = Nmax-N1.0
  N0 = rbind(N0, cbind(N1.0, N2.0) )
  data[, 1:2] = round(N0)
  return(data)
}

doe.grid = function(n,Nmax){
  data = matrix(NA, nrow=n, ncol=4)
  grid.axis = round((0:7)*Nmax/7)
  N0 = as.matrix(expand.grid(x=grid.axis, y=grid.axis)) # full grid
  N0[1, ] = c(grid.axis[2],grid.axis[2]) # replace origin
  data[, 1:2] = do.call("rbind", rep(list(N0), n/64))
  return(data)
}

doe.gridlog = function(n,Nmax){
  data = matrix(NA, nrow=n, ncol=4)
  grid.axis = round((1.5^(0:7)-1.5^0)*(Nmax)/(1.5^7-1.5^0))
  N0 = as.matrix(expand.grid(x=grid.axis, y=grid.axis)) # full grid
  N0[1, ] = c(grid.axis[2],grid.axis[2]) # replace origin
  data[, 1:2] = do.call("rbind", rep(list(N0), n/64))
  return(data)
}

doe.diags = function(n,Nmax){
  data = matrix(NA, nrow=n, ncol=4)
  data.diags = matrix(NA, nrow=1, ncol=2) # initialize matrix
  n.diags = floor(n/8)
  axis.diags = Nmax/8 * (1:8)
  for(i in 1:length(axis.diags)){
    diags = seq(from=0, to=axis.diags[i], length.out=n.diags)
    data.diags = rbind(data.diags, cbind(diags, rev(diags)))
  }
  data.diags = data.diags[-1, ] # delete initialization row
  data[, 1:2] = round(data.diags)
  return(data)
}

doe.random.3d = function(n,Nmax){
  data = matrix(NA, nrow=n, ncol=6) 
  for(i in 1:3){
    data[, i] = sample(0:Nmax, size=n, replace=TRUE)
  }
  return(data)
}

