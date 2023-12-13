

Holling2 = function(t, N, parms){
  with(parms,{ 
    dN = rep(0,m)
    denominator = 1.0
    for(i in 1:m) { denominator = denominator + (a[i]*h[i]*N[i]) }
    for(i in 1:m) { dN[i] = -a[i]*N[i] / denominator}
    return(list(dN))
  })
}

Holling3 = function(t, N, parms){
  with(parms,{ 
    dN = rep(0,m)
    denominator = 1.0
    for(i in 1:m) { denominator = denominator + (a[i]*h[i]*N[i]^(1+q)) }
    for(i in 1:m) { dN[i] = -a[i]*N[i]^(1+q) / denominator}
    return(list(dN))
  })
}

Yodzis = function(t, N, parms){
  with(parms,{ 
    dN = rep(0,m)
    denominator = 0.0
    for(i in 1:m){ denominator = denominator + w[i] * ( N[i]^r + a[i]*h[i]*N[i]^(1+r) ) }
    for(i in 1:m){ dN[i] = -a[i]*w[i]*N[i]^(1+r)/denominator }
    return(list(dN))
  })
}

Generalized = function(t, N, parms){
  with(parms,{
    dN = rep(0,m)
    denominator = sum(w[1:m]*N^r) + sum(a[1:m]*w[1:m]*h[1:m]*N^(1+q+r))
    for(i in 1:m) dN[i] = -a[i]*w[i]*N[i]^(1+q+r)/denominator
    return(list(dN))
  }) # end with
}

